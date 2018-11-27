#include "simulation.h"
#include <ctime>
#include <Eigen/Dense>
#include <random>
#include <iostream>
#include <random>
#include "triangulation.h"
#include <fstream>
#include "clock.h"
namespace tri
{
Simulation::Simulation(int n_cam, int n_pts, const double& sig)
    :n_cam_(n_cam), n_pts_(n_pts), err_sigma_(sig)
{
    K_ << 400, 0, 512,
            0, 400, 512,
            0, 0, 1;

    cam_h_ = 1024;
    cam_w_ = 1024;
    p_tri_ = nullptr;
    //std::cout << K_ << std::endl;
}

void Simulation::projection()
{
    assert(pts_.size() == n_pts_);
    assert(cams_.size() == n_cam_);

    pro_uvs_.resize(n_cam_);
    valid_proj_.resize(n_cam_);

    for(int i_cam = 0; i_cam < n_cam_; i_cam++)
    {
        //
        pro_uvs_[i_cam].resize(n_pts_);
        valid_proj_[i_cam].resize(n_pts_);

        const Camera & cur_cam = cams_[i_cam];


        for(int i_pt = 0; i_pt < n_pts_; i_pt ++)
        {
            bool flag = cur_cam.project(pts_[i_pt], pro_uvs_[i_cam][i_pt]);
            valid_proj_[i_cam][i_pt] = static_cast<int>(flag);
        }
    }
}

void Simulation::corrupt()
{
    //srand((unsigned)time(NULL));
    std::default_random_engine generator(time(NULL));;
    std::normal_distribution<double> distribution(0, err_sigma_);
    std::uniform_int_distribution<int> theta(-10000, 10000);

    corrupted_uvs_.resize(n_cam_);


    for(int i_cam = 0; i_cam < n_cam_; i_cam++)
    {
        corrupted_uvs_[i_cam].resize(n_pts_);

        for(int i_pt = 0; i_pt < n_pts_; i_pt ++)
        {
            if(valid_proj_[i_cam][i_pt])
            {
                double mag = distribution(generator);
                double angle =  M_PI * (theta(generator) / 10000.0);

                corrupted_uvs_[i_cam][i_pt] = pro_uvs_[i_cam][i_pt] + Eigen::Vector2d(mag * cos(angle), mag * sin(angle));
            }
        }
    }
}

void Simulation::genBearingVec()
{
    b_vec_.resize(n_cam_);

    for(int i_cam = 0; i_cam < n_cam_; i_cam++)
    {
        b_vec_[i_cam].resize(n_pts_);
        Camera& cur_cam = cams_[i_cam];

        for(int i_pt = 0; i_pt < n_pts_; i_pt ++)
        {
            if(valid_proj_[i_cam][i_pt])
            {
                Eigen::Vector2d uv_pro;
                b_vec_[i_cam][i_pt] = cur_cam.unproject(corrupted_uvs_[i_cam][i_pt]);
                cur_cam.project(b_vec_[i_cam][i_pt], uv_pro);
                //std::cout << b_vec_[i_cam][i_pt] << std::endl;
            }
        }
    }
    //std::cout << "bearing: " << b_vec_[n_cam_ - 1][n_pts_ - 1] << std::endl;
}

void Simulation::computeErr()
{
    assert(pts_r_.size() == n_pts_);

    n_valid_proj_ = std::vector<int>(n_pts_, 0);
    pt_err_vec_ = std::vector<double>(n_pts_, 0);

    reproj_err_.resize(n_cam_);


    for(int i_cam = 0; i_cam < n_cam_; i_cam++)
    {
        reproj_err_[i_cam].resize(n_pts_);

        Camera& cur_cam = cams_[i_cam];

        for(int i_pt = 0; i_pt < n_pts_; i_pt ++)
        {
            if(valid_proj_[i_cam][i_pt])
            {
                Eigen::Vector2d tmp_uv;
                bool flag = cur_cam.project(pts_r_[i_pt], tmp_uv);


                if(flag)
                {
                    tmp_uv -= corrupted_uvs_[i_cam][i_pt];

                    reproj_err_[i_cam][i_pt] = tmp_uv.norm();

                    n_valid_proj_[i_pt] ++;
                    pt_err_vec_[i_pt] += reproj_err_[i_cam][i_pt];

                } //if
                else
                {
                    //valid_proj_[i_cam][i_pt] = 0;
                }

            }  // if valid
        } // for i_pt
    } // for i_cam

    for(int i_pt = 0; i_pt < n_pts_; i_pt++)
    {
        if(n_valid_proj_[i_pt] > 0)
        {
            pt_err_vec_[i_pt] /= n_valid_proj_[i_pt];
        }
        else
        {
            n_valid_proj_[i_pt] = 0;
            pt_err_vec_[i_pt] = 0;
        }

    }

}

void Simulation::print_result()
{
    double sum_err = 0;
    double sum_pt = 0;
    for(int i = 0; i < n_pts_; i++)
    {
        if(n_valid_proj_[i] > 1)
        {
            sum_err += pt_err_vec_[i];
            sum_pt += 1;

        }
    }

    double mean_err = sum_err / sum_pt;

    std::cout << "Solver: " << p_tri_->name()
              << ", average reprojection error: " << mean_err
              << " px, runtime: " << run_time_ << "ms" << std::endl;
}

void Simulation::reconstruct(int alg)
{

    TriangulationMVMP init_alg;



    std::vector<Eigen::Vector3d> pts, b_vec;
    std::vector<Camera> cams;
    cams.reserve(n_cam_);
    pts.reserve(n_pts_);
    b_vec.reserve(n_pts_);
    pts_r_.resize(n_pts_);
    sucessed_.resize(n_pts_);

    if(alg == 0)
    {
        p_tri_ = std::make_shared<TriangulationMVMP>();
    }
    else if(alg == 1)
    {
        p_tri_ = std::make_shared<TriangulationIRMP>( );
    }
    else if(alg == 2)
    {
        p_tri_ = std::make_shared<TriangulationCeres>( );

    }



    Clock cl;
    cl.reset();

    run_time_ = 0;

    //static int run_pip = 1;
    // std::string file_name = "pip_test" + std::to_string(run_pip ++) + ".txt";

    // std::ofstream ofs(file_name);

    for(int i_pt = 0; i_pt < n_pts_; i_pt++)
    {
        pts.clear();
        cams.clear();
        b_vec.clear();
        //ofs << "pt: " << i_pt << std::endl;
        for(int i_cam = 0; i_cam < n_cam_; i_cam++)
        {
            if(valid_proj_[i_cam][i_pt])
            {
                cams.push_back(cams_[i_cam]);

                b_vec.push_back(b_vec_[i_cam][i_pt]);

             }
        }
        bool flag = false;
        if(cams.size() >= 2)
        {
            if(alg == 2)
            {

                init_alg.init(cams, b_vec);
                flag = init_alg.run(pts_r_[i_pt]);
                if(flag)
                {
                    p_tri_ -> init(cams, b_vec);



                    flag = p_tri_ -> run(pts_r_[i_pt]);
                }
            }
            else
            {
                p_tri_ -> init(cams, b_vec);
                flag = p_tri_ -> run(pts_r_[i_pt]);

            }
        }
        sucessed_[i_pt] = flag;
        //ofs << i_pt << ": " << flag << " " << pts_r_[i_pt] << std::endl;

    }
    run_time_ = cl.now() / 1000.0;
}

void SimulationConfig1::genSynData()
{
    // generate cams
    cams_.reserve(n_cam_);

    double dis = 10;
    double tiny_val = 0.003;
    double offset = - dis / 4;
    for(int i = 0; i < n_cam_; i++)
    {
        Eigen::Matrix3d R = Eigen::Matrix3d::Identity();
        Eigen::Vector3d t = Eigen::Vector3d::Zero();


        double scale = (n_cam_ - i) / 2.9;
        scale *= scale;

        t(0) = scale * tiny_val + offset;
        t(1) = scale * tiny_val + offset;
        t(2) = dis * i / n_cam_ ;

        t = -R.transpose() * t;
        R = R.transpose().eval();
        cams_.push_back(Camera(R, t, K_, cam_w_, cam_h_));
    }

    // generate pts;

    std::default_random_engine generator(time(NULL));
    const double scale = 1000;
    std::uniform_int_distribution<int> dist_xy(- scale * dis / 2, scale * dis / 2);

    std::uniform_int_distribution<int> dist_z(scale * dis * 1.1, scale * dis * 1.15);

    pts_.reserve(n_pts_);
    for (int i = 0; i < n_pts_; i++)
    {
        pts_.push_back(Eigen::Vector3d(dist_xy(generator) / scale, dist_xy(generator) / scale, dist_z(generator) / scale));
    }
}


void SimulationConfig2::genSynData()
{
    // generate cams
    cams_.reserve(n_cam_);

    double dis = 10;
    double tiny_val = 0.03;
    double offset = - dis / 4;
    for(int i = 0; i < n_cam_; i++)
    {
        Eigen::Matrix3d R = Eigen::Matrix3d::Identity();
        Eigen::Vector3d t = Eigen::Vector3d::Zero();


        double scale = (n_cam_ - i) / 2.9;
        scale *= scale;

        t(0) = sin(i * tiny_val) * dis / 3;
        t(1) = cos(i * tiny_val) * dis / 3;
        t(2) = dis * i / n_cam_ ;

        t = -R.transpose() * t;
        R = R.transpose().eval();
        cams_.push_back(Camera(R, t, K_, cam_w_, cam_h_));
    }

    // generate pts;

    std::default_random_engine generator(time(NULL));
    const double scale = 1000;
    std::uniform_int_distribution<int> dist_xy(- scale * dis / 2, scale * dis / 2);

    std::uniform_int_distribution<int> dist_z(- scale * dis * 0.1, scale * dis * 1.1);

    pts_.reserve(n_pts_);
    for (int i = 0; i < n_pts_; i++)
    {
        pts_.push_back(Eigen::Vector3d(dist_xy(generator) / scale, dist_xy(generator) / scale, dist_z(generator) / scale));
    }
}


void SimulationConfig3::genSynData()
{
    // generate cams
    cams_.reserve(n_cam_);

    double dis = 10;
    double tiny_val = 0.003;
    double offset = - dis / 4;

    double radius = 10;

    for(int i = 0; i < n_cam_; i++)
    {

        double theta = 2 * M_PI * i / n_cam_;

        Eigen::Vector3d z_new(cos(theta), sin(theta), 0);

        double angle_1 = - 45 * M_PI / 180;

        Eigen::Quaterniond q1(cos(angle_1), 0, sin(angle_1), 0);

        double angle_2 = theta / 2 + M_PI / 2;

        Eigen::Quaterniond q2(cos(angle_2), sin(angle_2), 0, 0);


        Eigen::Matrix3d R = q2.toRotationMatrix() * q1.toRotationMatrix();


        Eigen::Vector3d t = z_new * radius;

        t = -R.transpose() * t;
        R = R.transpose().eval();
        cams_.push_back(Camera(R, t, K_, cam_w_, cam_h_));
    }

    // generate pts;

    std::default_random_engine generator(time(NULL));
    const double scale = 1000;
    double r = radius * 0.3;
    std::uniform_int_distribution<int> dist(- scale * r, scale * r);

    pts_.reserve(n_pts_);
    for (int i = 0; i < n_pts_; i++)
    {
        pts_.push_back(Eigen::Vector3d(dist(generator) / scale, dist(generator) / scale, dist(generator) / scale));
    }
}


void SimulationConfig4::genSynData()
{

    // generate cams
    cams_.reserve(n_cam_);
    std::default_random_engine generator(time(NULL));
    const double scale = 1000;
    std::uniform_int_distribution<int> ran_uniform_4(- scale * 4 , scale * 4 );
    std::uniform_int_distribution<int> ran_uniform_5(- scale * 5 , scale * 5 );

    double dis = 10;
    double tiny_val = 0.003;
    double offset = - dis / 4;
    for(int i = 0; i < n_cam_; i++)
    {
        Eigen::Vector4d quat_vec(ran_uniform_4(generator) / scale,ran_uniform_4(generator) / scale,ran_uniform_4(generator) / scale,ran_uniform_4(generator) / scale);
        quat_vec.normalize();
        Eigen::Quaterniond quat(quat_vec);

        Eigen::Matrix3d R = quat.toRotationMatrix();
        Eigen::Vector3d t (ran_uniform_4(generator) / scale, ran_uniform_4(generator) / scale, ran_uniform_4(generator) / scale);

        t = -R.transpose() * t;
        R = R.transpose().eval();
        cams_.push_back(Camera(R, t, K_, cam_w_, cam_h_));
    }

    // gen pts
    pts_.reserve(n_pts_);
    for (int i = 0; i < n_pts_; i++)
    {
        pts_.push_back(Eigen::Vector3d(ran_uniform_5(generator) / scale, ran_uniform_5(generator) / scale, ran_uniform_5(generator) / scale));
    }
}
}

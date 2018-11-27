#include "triangulation.h"
#include <Eigen/Cholesky>
#include <Eigen/Eigenvalues>
#include <Eigen/Jacobi>
#include <fstream>
#include <ceres/ceres.h>
namespace tri
{
void TriangulationMVMP::init(const std::vector<Camera> &cams,
                             const std::vector<Eigen::Vector3d> &b)
{
    cams_ = cams;
    b_vec_ = b;
    name_ = "MVMP";
}

bool TriangulationMVMP::run(Eigen::Vector3d& res)
{
    //std::vector<Eigen::Matrix3d> A_vec;
    //A_vec.reserve(cams_.size());
    Eigen::Matrix3d A_sum = Eigen::Matrix3d::Zero();
    Eigen::Vector3d b_sum = Eigen::Vector3d::Zero();
    for(int i = 0; i < b_vec_.size(); i ++)
    {
        Eigen::Vector3d cur_b = b_vec_[i];
        Eigen::RowVector3d v_t = cur_b.transpose();

      //  std::cout << cur_b << std::endl;

        Eigen::Matrix3d A = Eigen::Matrix3d::Identity() - cur_b * v_t;
      //  std::cout << A << std::endl;

        A_sum += A;
        b_sum += (A * cams_[i].t_inv());

    }

    res = A_sum.jacobiSvd(Eigen::ComputeFullU | Eigen::ComputeFullV).solve(b_sum);


    if(res.allFinite())
        return true;
    else
        return false;
}


void TriangulationIRMP::init(const std::vector<Camera> &cams,
                             const std::vector<Eigen::Vector3d> &b)
{
    cams_ = cams;
    b_vec_ = b;

    name_ = "IRMP";
}




bool TriangulationIRMP::run(Eigen::Vector3d& res)
{

    std::vector<Eigen::Matrix3d> A_vec(b_vec_.size());
    std::vector<Eigen::Vector3d> t_vec(b_vec_.size());
    for(int i = 0; i < b_vec_.size(); i ++)
    {
        Eigen::Vector3d cur_b = b_vec_[i];
        Eigen::RowVector3d v_t = cur_b.transpose();

        A_vec[i] = Eigen::Matrix3d::Identity() - cur_b * v_t;
        t_vec[i] = A_vec[i] * cams_[i].t_inv();
    }

    // init by mvmp
    Eigen::Matrix3d A_sum = Eigen::Matrix3d::Zero();
    Eigen::Vector3d b_sum = Eigen::Vector3d::Zero();
    for(int i = 0; i < b_vec_.size(); i ++)
    {
        Eigen::Matrix3d& A = A_vec[i];
        A_sum += A;
        b_sum += t_vec[i];
    }

    res = A_sum.jacobiSvd(Eigen::ComputeFullU | Eigen::ComputeFullV).solve(b_sum);
    if(!res.allFinite())
    {
        return false;
    }

    Eigen::Vector3d pre_pt;
    Eigen::Vector3d &pt = res;


    for(int iter = 1; iter < max_iters_; iter++)
    {
        pre_pt = pt;

        A_sum.setZero();
        b_sum.setZero();


        for(int i = 0; i < b_vec_.size(); i ++)
        {
            Eigen::Vector3d op = pre_pt - cams_[i].t_inv();
            double cur_w = 1 / op.dot(op);
            Eigen::Matrix3d& A = A_vec[i];
            Eigen::Vector3d sine_op = A * op;
            double cur_err = cur_w * sine_op.dot(sine_op);


            A_sum += A * cur_w;
            b_sum += cur_w * (t_vec[i] + cur_err * op );
        }


        pt = A_sum.jacobiSvd(Eigen::ComputeFullU | Eigen::ComputeFullV).solve(b_sum);
        if(!pt.allFinite())
        {
            return false;
        }
        if((pt - pre_pt).norm() < 0.001)
            break;
    }

    return true;

}


void TriangulationCeres::init(const std::vector<Camera> &cams, const std::vector<Eigen::Vector3d> &uvs)
{
    cams_ = cams;

    p_pro_ = new ceres::Problem;

    uv_vec_.resize(uvs.size());
    for(int i = 0; i < uvs.size(); i ++)
    {
        Eigen::Vector3d uv_c = cams[i].R() * uvs[i];

        uv_vec_[i](0) = uv_c(0) / uv_c(2);
        uv_vec_[i](1) = uv_c(1) / uv_c(2);
        ProjectionFactor * factor = new ProjectionFactor(uv_vec_[i], cams_[i]);
        p_pro_->AddResidualBlock(factor, NULL, state_);

    }


    name_ = "Ceres";
}

bool TriangulationCeres::run(Eigen::Vector3d &res)
{


    state_[0] = res(0);
    state_[1] = res(1);
    state_[2] = res(2);

    ceres::Solver::Options options;

    //options.linear_solver_type = ceres::DENSE_SCHUR;
    //options.trust_region_strategy_type = ceres::LEVENBERG_MARQUARDT;
    options.max_num_iterations = 10;
    //options.minimizer_progress_to_stdout = true;
    ceres::Solver::Summary summary;
    ceres::Solve(options, p_pro_, &summary);


    res(0) = state_[0];
    res(1) = state_[1];
    res(2) = state_[2];

    //ceres
}


int ProjectionFactor::id_count_ = 0;

ProjectionFactor::ProjectionFactor(const Eigen::Vector2d& uv, const Camera& cam )
    :uv_(uv), cam_(cam), id_(id_count_++)
{
    ;
}

bool ProjectionFactor::Evaluate(double const *const *parameters, double *residuals, double **jacobians) const
{
    Eigen::Vector3d pt(parameters[0][0], parameters[0][1], parameters[0][2]);

    Eigen::Vector3d pt_c = cam_.w2c(pt);

    Eigen::Vector2d uv = pt_c.head(2) / pt_c(2);

    Eigen::Map<Eigen::Vector2d> residual(residuals);
    residual = uv - uv_;

  // std::cout << id_ << " " << residual << std::endl;


    // compute jacobian
    if (jacobians)
    {
        Eigen::Matrix3d R = cam_.R();


        double z_inv = 1.0 / pt_c(2);
        double z_inv2 = z_inv * z_inv;

        double x_z2 = pt_c(0) * z_inv2;
        double y_z2 = pt_c(1) * z_inv2;


        if (jacobians[0])
        {
            Eigen::Map<Eigen::Matrix<double, 2, 3, Eigen::RowMajor>> jac_mat(jacobians[0]);


            jac_mat(0, 0) = z_inv;
            jac_mat(0, 1) = 0.0;
            jac_mat(0, 2) = -x_z2;
            jac_mat(1, 0) = 0.0;
            jac_mat(1, 1) = z_inv;
            jac_mat(1, 2) = -y_z2;

            jac_mat = jac_mat * R;

            //std::cout << jac_mat << std::endl;
        }


    }

    return true;
}

}

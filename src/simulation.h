#pragma once
#include "camera.h"
#include <vector>
#include "triangulation.h"


namespace tri
{
class Simulation
{
public:
    explicit Simulation(int n_cam = 100, int n_pts = 5000, const double& sig = 10);

    //! Generate synthetic data.
    virtual void genSynData() = 0;

    //! Project 3D pt to the camera plane.
    void projection();

    //! corrupt
    void corrupt();

    //! Generate bearing vector.
    void genBearingVec();

    //!
    void reconstruct(int alg);


    void computeErr();

    void print_result();

protected:

    //!
    int n_cam_;
    int n_pts_;


    //! Synthetic cams
    std::vector<Camera> cams_;

    //! Synthetic pts
    std::vector<Eigen::Vector3d> pts_;

    //! Projection result size = n_cam_ * n_pts_
    std::vector<std::vector<Eigen::Vector2d> > pro_uvs_;

    //! If this is a valid projection
    std::vector<std::vector<int> > valid_proj_;

    //! Corrupted uvs.
    std::vector<std::vector<Eigen::Vector2d> > corrupted_uvs_;

    //! Bearing vector.
    std::vector<std::vector<Eigen::Vector3d> > b_vec_;



    //! Reconstructed pts
    std::vector<Eigen::Vector3d> pts_r_;
    std::vector<int> sucessed_;


    //! reprojection errors
    std::vector<std::vector<double> > reproj_err_;
    //! number of valid projection
    std::vector<int> n_valid_proj_;
    std::vector<double> pt_err_vec_;


    Eigen::Matrix3d K_;
    int cam_w_;
    int cam_h_;
    double run_time_;

    double err_sigma_;
    TriangulationPtr p_tri_;

};


class SimulationConfig1 : public Simulation
{
public:
    explicit SimulationConfig1(int n_cam = 100, int n_pts = 5000, const double& sig = 10) : Simulation(n_cam, n_pts, sig){;}
    virtual void genSynData();

private:
};

class SimulationConfig2 : public Simulation
{
public:
    explicit SimulationConfig2(int n_cam = 100, int n_pts = 5000, const double& sig = 10) : Simulation(n_cam, n_pts, sig){;}
    virtual void genSynData();

private:
};


class SimulationConfig3 : public Simulation
{
public:
    explicit SimulationConfig3(int n_cam = 100, int n_pts = 5000, const double& sig = 10) : Simulation(n_cam, n_pts, sig){;}
    virtual void genSynData();

private:
};


class SimulationConfig4 : public Simulation
{
public:
    explicit SimulationConfig4(int n_cam = 100, int n_pts = 5000, const double& sig = 10) : Simulation(n_cam, n_pts, sig){;}
    virtual void genSynData();
};
}

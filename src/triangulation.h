#pragma once
#include <Eigen/Core>
#include "camera.h"
#include "ceres/ceres.h"
#include <memory>
namespace tri
{
class Triangulation
{
public:
    Triangulation():max_iters_(10){}
    virtual void init(const std::vector<Camera>& cams,
              const std::vector<Eigen::Vector3d>& b
              ) = 0;

    virtual bool run(Eigen::Vector3d& res) = 0;

    inline std::string name() const {return name_;}
protected:
    int max_iters_;
    std::string name_;
};

typedef std::shared_ptr<Triangulation> TriangulationPtr;


class TriangulationMVMP : public Triangulation
{
public:
    void init(const std::vector<Camera>& cams,
                      const std::vector<Eigen::Vector3d>& pts
                      ) final;
    bool run(Eigen::Vector3d& res) final;

private:

    std::vector<Camera> cams_;
    std::vector<Eigen::Vector3d> b_vec_;

};

class TriangulationIRMP : public Triangulation
{
public:
    void init(const std::vector<Camera>& cams,
                      const std::vector<Eigen::Vector3d>& pts
                      ) final;
    bool run(Eigen::Vector3d& res) final;

private:

    std::vector<Camera> cams_;
    std::vector<Eigen::Vector3d> b_vec_;
};


class TriangulationNN : public Triangulation
{
public:
    void init(const std::vector<Camera>& cams,
                      const std::vector<Eigen::Vector3d>& pts
                      ) final;
    bool run(Eigen::Vector3d& res) final;

private:

    std::vector<Camera> cams_;
    std::vector<Eigen::Vector3d> b_vec_;
};


class TriangulationCeres : public Triangulation
{
public:
    void init(const std::vector<Camera>& cams,
                      const std::vector<Eigen::Vector3d>& uvs
                      ) final;
    bool run(Eigen::Vector3d& res) final;

    ~TriangulationCeres(){if(!p_pro_) delete p_pro_;}
private:


    ceres::Problem* p_pro_;

    std::vector<Camera> cams_;
    std::vector<Eigen::Vector2d> uv_vec_;
    double state_[3];
};


class ProjectionFactor : public ceres::SizedCostFunction<2, 3>
{
public:
    ProjectionFactor(const Eigen::Vector2d& uv, const Camera& cam );
    virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const;

private:
    Eigen::Vector2d uv_;
    Camera cam_;
    int id_;
    static int id_count_;
};

}

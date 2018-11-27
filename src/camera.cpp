#include "camera.h"
namespace tri
{
Camera::Camera(const Eigen::Matrix3d& att,
               const Eigen::Vector3d& loc,
               const Eigen::Matrix3d& K,
               int w,
               int h)
    :R_(att), t_(loc), K_(K), w_(w), h_(h)
{
    K_inv_ = K_.inverse();
    t_inv_ = - R_.transpose() * t_;
}


Eigen::Vector3d Camera::unproject(const Eigen::Vector2d &uv) const
{
    Eigen::Vector3d pt{0, 0, 1};
    pt(0) = uv(0);
    pt(1) = uv(1);
    pt(2) = 1;
    pt = R_.transpose() * K_inv_ * pt;
    return pt / pt.norm();
}
Eigen::Vector3d Camera::w2c(const Eigen::Vector3d& pt_w) const
{
    return R_ * pt_w + t_;
    //Eigen::Vector3d pt_c = ;

    //return pt_c.head(2) / pt_c(2);
}

bool Camera::project(const Eigen::Vector3d &pt_w, Eigen::Vector2d &uv) const
{
    Eigen::Vector3d pt_c = R_ * pt_w + t_;
    if(pt_c(2) < 0)
        return false;

    pt_c = K_ * pt_c;

    uv = pt_c.head(2) / pt_c(2);

    if(uv(0) < 0 || uv(0) > w_ - 1 ||
            uv(1) < 0 || uv(1) > h_ - 1)
    {
        return false;
    }

    return true;
}
}

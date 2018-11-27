#pragma once

#include <Eigen/Core>
#include <Eigen/LU>
//#include <Eigen/Geometry>

namespace tri
{
    class Camera
    {
    public:

        explicit Camera(const Eigen::Matrix3d& att,
                        const Eigen::Vector3d& loc,
                        const Eigen::Matrix3d& K,
                        int w,
                        int h);


        //! project pt (world coordinate) to the camera plane
        //! return false if failed.
        bool project(const Eigen::Vector3d& pt_w, Eigen::Vector2d& uv) const;

        Eigen::Vector3d w2c(const Eigen::Vector3d& pt_w) const;



        Eigen::Vector3d unproject(const Eigen::Vector2d& uv) const;


        inline Eigen::Matrix3d R() const {return R_;}
        //inline Eigen::Matrix3d kR() const {return K_ * R_;}
        inline Eigen::Vector3d t() const {return t_;}
        inline Eigen::Vector3d t_inv() const {return t_inv_;}
    private:
        //! attitude
        Eigen::Matrix3d R_;

        //! location
        Eigen::Vector3d t_;

        Eigen::Vector3d t_inv_;

        //! inner parameters.
        Eigen::Matrix3d K_;

        //! inverse of K_
        Eigen::Matrix3d K_inv_;

        //! width
        int w_;

        //! height
        int h_;

    };
}

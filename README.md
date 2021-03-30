# Large-Scale Trajectory Optimizer
This is __probably the fastest__ minimum jerk or minimum snap trajectory generation you can find.

It also provides __analytical gradient__ of energy with respect to time allocations and waypoints.

## 1. About
This is a highly-optimized implementation for minimum jerk/snap trajectories with exact gradient w.r.t. time allocation and waypoints. It is based on completely analytical results of the following paper.

Related Paper: [__Generating Large-Scale Trajectories with Polynomial Double Descriptions__](https://zhepeiwang.github.io/pubs/icra_2021_largescale.pdf)

(Accepted by ICRA2021. [Video](https://www.youtube.com/watch?v=tA3fIyggH4I))

If this repo helps you, please cite our paper:

    @article{WANG2020GLST,
        title={Generating Large-Scale Trajectories Efficiently using Double Descriptions of Polynomials},
        author={Wang, Zhepei and Ye, Hongkai and Xu, Chao and Gao, Fei},
        journal={arXiv preprint arXiv:2011.02662v2},
        year={2020}
    }

__Author__: [Zhepei Wang](https://zhepeiwang.github.io/) and [Fei Gao](https://ustfei.com/) from the [ZJU Fast Lab](http://zju-fast.com/).

## 2. Examples

__Example 1__ gives the computation speed of our implementation.

__Example 2__ will be released soon. It is a high-performance large-scale trajectory optimizer. It achieves almost the same trajectory quality as the global trajectory optimizer in [Teach-Repeat-Replan](https://github.com/HKUST-Aerial-Robotics/Teach-Repeat-Replan) while using significantly less computation time.

## 3. Performance

We compare our implementation with four existing works. The performance is shown as follows.

<p align="center">
  <img src="misc/ModerateScale.png" width = "640" height = "375"/>
</p>
<p align="center">
  <img src="misc/LargeScale.png" width = "640" height = "272"/>
</p>

<sub> <em> Burke et al., “Generating minimum-snap quadrotor trajectories really fast,” IROS 2020. </em> </sub> <br/>
<sub> <em> Bry et al., “Aggressive flight of fixed-wing and quadrotor aircraft in dense indoor environments,” IJRR 2015. </em> </sub> <br/>
<sub> <em> Mellinger et al., “Minimum snap trajectory generation and control for quadrotors,” ICRA 2011. </em> </sub>

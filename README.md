# Large-Scale Trajectory Optimizer
This is probably __the fastest__ minimum jerk or minimum snap trajectory generation you can find.

It also provides __analytical gradient__ for jerk/snap energy w.r.t. time allocations and waypoints.

## 0. About
If this repo helps you, please cite the paper below.

Related Paper: __Generating Large-Scale Trajectories with Polynomial Double Descriptions__

(Accepted by ICRA2021. [ArXiv version](https://arxiv.org/abs/2011.02662v2)and [video](https://zhepeiwang.github.io/pubs/icra_2021_sub_largescale.mp4) are all available now.)

__Author__: [Zhepei Wang](https://zhepeiwang.github.io/) and [Fei Gao](https://ustfei.com/) from the [ZJU Fast Lab](http://zju-fast.com/).

## 1. Examples

__Example 1__ gives the computation speed of our implementation.

__Example 2__ will be released soon. It is a high-performance large-scale trajectory optimizer. It achieves almost the same trajectory quality as the global trajectory optimizer in [Teach-Repeat-Replan](https://github.com/HKUST-Aerial-Robotics/Teach-Repeat-Replan) while using significantly less computation time.

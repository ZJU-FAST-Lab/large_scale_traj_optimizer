# Large-Scale Trajectory Optimizer
This is __probably the fastest__ minimum jerk or minimum snap trajectory generation you can find.

It also provides __analytical gradient__ of energy with respect to time allocations and waypoints.

## 0. About
If this repo helps you, please cite the paper below.
This is a highly-optimized implementation for minimum jerk/snap trajectories with exact gradient w.r.t. time allocation and waypoints. It is based on completely analytical results of the following paper. If this repo helps you, please cite our paper:

Related Paper: __Generating Large-Scale Trajectories with Polynomial Double Descriptions__

(Accepted by ICRA2021. [ArXiv version](https://arxiv.org/abs/2011.02662v2) and [video](https://www.youtube.com/watch?v=tA3fIyggH4I) are all available now.)

__Author__: [Zhepei Wang](https://zhepeiwang.github.io/) and [Fei Gao](https://ustfei.com/) from the [ZJU Fast Lab](http://zju-fast.com/).

## 1. Examples

__Example 1__ gives the computation speed of our implementation.

__Example 2__ will be released soon. It is a high-performance large-scale trajectory optimizer. It achieves almost the same trajectory quality as the global trajectory optimizer in [Teach-Repeat-Replan](https://github.com/HKUST-Aerial-Robotics/Teach-Repeat-Replan) while using significantly less computation time.

# LBFGS-Lite

A header-only LBFGS unconstrained optimizer.

## 笔记

```C++
inline int lbfgs_optimize(int n,
        double *x,
        double *ptr_fx,
        lbfgs_evaluate_t proc_evaluate,
        lbfgs_stepbound_t proc_stepbound,
        lbfgs_progress_t proc_progress,
        void *instance,
        lbfgs_parameter_t *_param)
```

调用函数`lbfgs_optimize`来进行优化，其参数包括：

- `n`：变量的个数
- `x`：变量数组。客户端可以根据该数组设定优化的初始变量
- `ptr_fx`：根据变量返回优化函数的最终结果，如果不需要优化函数，则可以设定为`NULL`
- `proc_evaluate`：在给定一组当前变量值的时候，该回调函数可以评估函数的梯度值，客户端必须实现与`lbfgs_evaluate_t`兼容的回调函数，并且将该函数指针传给该变量
- `proc_stepbound`：该回调函数提供搜索变量步长的上限值，提供线性搜索前变量的起始值，以及当前步长向量（可能是负梯度）。客户端可以实现此功能以实现更高效的搜索，如果不需要设置为`NULL`
- `proc_progress`：该回调函数接收最小化过程的进度（迭代次数，目标函数的当前值）。如果不需要设置进度报告，可以设置为`NULL`
- `instance`：客户端的实例指针，回调函数将接受此参数的值
- `_param`：指向表示 LBFGS 优化参数的结构的指针。客户端可以使用 NULL 以使用默认参数，可以通过调用函数`lbfgs_load_default_parameters()`来初始化默认参数的结构体
- 返回值：返回状态，如果为 0 则无任何错误，否则出错

```C++
struct lbfgs_parameter_t
{
    int mem_size;
    double g_epsilon;
    int past;
    double delta;
    int max_iterations;
    int max_linesearch;
    double min_step;
    double max_step;
    double f_dec_coeff;
    double s_curv_coeff;
    int abs_curv_cond;
    double xtol;
    int shrink_type;
}
```

优化器的参数如上结构体所示，其中参数的含义为：

- `mem_size`：近似逆 hessian 矩阵的修正次数。L-BFGS 例程存储前 m 次迭代的计算结果以近似当前迭代的逆 Hessian 矩阵。此参数控制有限内存（修正）的大小，默认值为 `8` ，不建议使用小于 3 的值，较大的值会导致计算时间过长。
- `g_epsilon`：用于梯度收敛测试的 Epsilon 。该参数决定了求解精度，当满足条件 `||g|| < g_epsilon * max(1, ||x||)` 时最小化终止。其中 `||.||` 表示欧几里得范数，默认值为 `1e-5`
- `past`：基于增量的收敛测试的距离。此参数确定迭代中的距离，以计算目标函数的下降率。如果此参数的值为零，则库不执行基于增量的收敛测试。默认值为 `0`。
- `delta`：用于收敛测试的 Delta。该参数决定了目标函数的最小下降率。当满足以下条件时，库停止迭代：`|f' - f| / f < delta`，其中 `f'` 是过去迭代之前的目标值，`f` 是当前迭代的目标值。默认值为 `1e-5`。
- `max_iterations`：最大迭代次数。当迭代计数超过此参数时，`lbfgs_optimize()` 函数使用 `LBFGSERR_MAXIMUMITERATION` 状态代码终止优化过程。 将此参数设置为零会继续优化过程，直到收敛或出错。 默认值为 `0`。
- `max_linesearch`：线搜索的最大试验次数。此参数控制线搜索例程每次迭代的函数和梯度评估的数量。默认值为 `60`。
- `min_step`：线搜索例程的最小步骤。默认值为 `1e-20`。除非指数对于正在使用的机器来说太大，或者除非问题非常严重（在这种情况下应该增加指数），否则不需要修改该值。
- `max_step`：线搜索的最大步长。默认值为 `1e+20`。除非指数对于正在使用的机器来说太大，或者除非问题非常严重（在这种情况下应该增加指数），否则不需要修改该值。
- `f_dec_coeff`：控制线搜索例程精度的参数。默认值为 `1e-4`。此参数应大于 `0` 且小于 `0.5`。
- `s_curv_coeff`：控制线搜索例程精度的参数。默认值为 `0.9`。如果函数和梯度评估相对于迭代的成本来说并不昂贵（有时在解决非常大的问题时就是这种情况），将这个参数设置为一个小的值可能是有利的。典型的小值是 `0.1`。该参数应大于 `f_dec_coeff` 参数 `(1e-4)` 且小于 `1.0`。
- `abs_curv_cond`：用于确定要使用的曲率条件的参数。默认值为 `1`，表示强 Wolfe 条件。如果值为 `0`，则使用弱 Wolfe 条件。
- `xtol`：浮点值的机器精度。默认值为 `1e-16`。此参数必须是客户端程序设置的正值以估计机器精度。如果不确定区间的相对宽度小于此参数，则行搜索例程将以状态代码 `(LBFGSERR_ROUNDING_ERROR)` 终止。
- `shrink_type`：在线搜索中确定区间收缩策略的参数。默认值为 `0`，保证随着试验次数 `m` 的增加，不确定性至少减少 `2^(-m/2)`。 More 和 Thuente 使用它。 如果值为 `1`，则保证至少 `2^(-m)` 降低不确定性。 Hager 和 Zhang 在 TOMS851 中使用它。 前者适合条件良好的函数，后者适合僵硬的函数。

## 0. About

**LBFGS-Lite** is a **C/C++** [**header-only**](https://en.wikipedia.org/wiki/Header-only) library for **unconstrained optimization** on **twice continuously differentiable functions**. The code is modified from [**liblbfgs**](https://github.com/chokkan/liblbfgs). Only necessary part is preserved for simplicity. Some engineering considerations are also added for improved robustness.

## 1. Features

- Only one header file "lbfgs.hpp" is required for usage.

- No dependencies except C/C++ standard library.

- The library is an implementation of [**Limited-Memory Broyden-Fletcher-Goldfarb-Shanno**](https://doi.org/10.1007/BF01589116) (LBFGS) with [**More-Thuente Line Search**](https://doi.org/10.1145/192115.192132) to ensure linear time/space complexity and [strong Wolfe conditions](https://en.wikipedia.org/wiki/Wolfe_conditions) in each iteration.

- The objective function is required to be **twice continuously differentiable** on its domain.

- The library provides an additional callback to utilize externally provided maximum feasible stepsize. It can be helpful when the [function is closed on a bounded open domain](https://en.wikipedia.org/wiki/Closed_convex_function) instead of the whole Euclidean space. The callback avoids function evaluations at infeasible region. This can help a lot for closed functions as long as the Newton step is not always clipped.

- Engineering features such as skipping update at extremely small curvature and lower safeguarding for ill-conditioned cases are also adopted.

- Instruction set dependent parts and L1 regularization parts in original code are removed. Multiple files are reorgainzed here with some additional modification.

## 2. How to use

See "lbfgs_sample.cpp" for details.

## 3. Planned features

- LBFGS is only proved to be globally convergent (convergent to stationaries for any initial guess) under convexity assumption, while it also works in many nonconvex cases and yields almost the best results. Slight modification in updating can ensure global convergence without convexity assumption, following a work by [Fukushima and Li](<https://doi.org/10.1016/S0377-0427(00)00540-9>).

- Although More-Thuente line search is already good enough, we plan to compare it with [Hager-Zhang line search](https://doi.org/10.1137/030601880) and [Nocedal's zoom line search](https://link.springer.com/book/10.1007%2F978-0-387-40065-5). Hager-Zhang linear search is reported to be full machine accuracy while More-Thuente's method can only attain the half. Nocedal's zoom is also reported to be reliable. We will compare them and choose the best one if the improvement is good enough.

- Further code optimization (instruction set independent) and reorganization.

## 6. Licence

LBFGS-Lite is distributed under the term of the MIT license according the previous version by Okazaki and the initial FORTRAN version by Nocedal. Please refer to LICENSE file in the distribution.

## 7. Maintaince

If any bug, please contact [Zhepei Wang](https://zhepeiwang.github.io/) (<wangzhepei@live.com>).

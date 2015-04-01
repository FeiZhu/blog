---
layout: post
title: Solve for Harmonic Displacements with FEM Discretization
---

For comparison reasons, I need to implement the SIGGRAPH2009 paper [Numerical Coarsening of Inhomogeneous Elastic Materials](http://users.cms.caltech.edu/~owhadi/publications/KMOD09.pdf). Everything else was set up except the harmonic displacements $ h_{\alpha\beta} (1 \le \alpha \le \beta \le d) $:

$$
f(z) =
 \left\{ 
 \begin{array}{rcl}
 \overline{\overline{z^2}+\cos z} & \mbox{for}
 & |z|<3 \\ 0 & \mbox{for} & 3\leq|z|\leq5 \\
 \sin\overline{z} & \mbox{for} & |z|>5
 \end{array}\right.
$$

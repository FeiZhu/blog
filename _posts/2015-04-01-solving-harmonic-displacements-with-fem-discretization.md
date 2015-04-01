---
layout: post
title: Solve for Harmonic Displacements with FEM Discretization
---

For comparison reasons, I need to implement the SIGGRAPH2009 paper [Numerical Coarsening of Inhomogeneous Elastic Materials](http://users.cms.caltech.edu/~owhadi/publications/KMOD09.pdf). Everything else was set up except the harmonic displacements $ h_{\alpha\beta} (1 \le \alpha \le \beta \le d) $:

$$
f(z) =
\begin{array}
div(\mathbf{C}:\epsilon(h\_{\alpha\beta})=0 \quad & inside \Omega \\
(\mathbf{C}:\epsilon(h\_{\alpha\beta}) \dot \mathbf{n} = \epsilon (x\_{\alpha}\mathbf{e}_{\beta})\dot \mathbf{n} \quad & for \mathbf{x} \in \partial \Omega
\end{array}
$$

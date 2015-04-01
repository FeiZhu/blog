---
layout: post
title: Solve for Harmonic Displacements with FEM Discretization
---

For comparison reasons, I need to implement the SIGGRAPH2009 paper [Numerical Coarsening of Inhomogeneous Elastic Materials](http://users.cms.caltech.edu/~owhadi/publications/KMOD09.pdf). Everything else was set up except the harmonic displacements $ \mathbf{h}_{\alpha\beta} (1 \le \alpha \le \beta \le d) $:

$$
\begin{cases}
&\textrm{div}(\mathbf{C}:\epsilon(\mathbf{h}_{\alpha\beta}))=0 \quad & \textrm{inside}\  \Omega \\
&(\mathbf{C}:\epsilon(\mathbf{h}_{\alpha\beta}))\cdot\mathbf{n}=\epsilon(x_{\alpha}\mathbf{e}_{\beta})\cdot\mathbf{n} \quad & \textrm{for} \ \mathbf{x} \in \partial \Omega
\end{cases}
$$

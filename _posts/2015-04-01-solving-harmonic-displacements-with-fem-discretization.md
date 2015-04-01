---
layout: post
title: Solve for Harmonic Displacements with FEM Discretization
use_math: true
---

For comparison reasons, I need to implement the SIGGRAPH2009 paper [Numerical Coarsening of Inhomogeneous Elastic Materials](http://users.cms.caltech.edu/~owhadi/publications/KMOD09.pdf). Everything else was set up except the harmonic displacements $ \mathbf{h}_{\alpha\beta} (1 \le \alpha \le \beta \le d) $:

$$
\begin{cases}
\textrm{div}(\mathbf{C}:\epsilon(\mathbf{h}_{\alpha\beta}))=0 \quad  &\textrm{inside}\  \Omega \\
(\mathbf{C}:\epsilon(\mathbf{h}\_{\alpha\beta}))\cdot\mathbf{n}=\epsilon(x\_{\alpha}\mathbf{e}_{\beta})\cdot\mathbf{n} \quad    &\textrm{for} \ \mathbf{x} \in \partial \Omega
\end{cases}
$$

This is in essence a static equivalence problem with Neumann boundary conditions prescribing surface tractions equal to $ \epsilon(x_{\alpha}\mathbf{e}_{\beta})\cdot\mathbf{n} $. Here we consider solving this equation with FEM discretization.

**Disclaimer:**  this article has loose terminology .
**Notation:** [Einstein summation convention](http://en.wikipedia.org/wiki/Einstein_notation) is used.

First we construct the weak form of the original problem  by multiplying both sides with the variation of displacement $ \delta\mathbf{u}$ and integrating over the domain:

$$
\int_{\Omega}\sigma_{ij,j}\delta\mathbf{u}_id\Omega = 0
$$

From tensor calculus we know that $\sigma\_{ij}\delta\mathbf{u}\_j = \sigma\_{ij}\delta\mathbf{u}\_{i,j}+\sigma\_{ij,j}\delta\mathbf{u}\_i $. Substitute into the equation and apply the [divergence theorem](http://en.wikipedia.org/wiki/Divergence_theorem), we get:

$$
\int_\Omega\sigma_{ij,j}\delta\mathbf{u}_id\Omega=-\int_\Omega\sigma_{ij}\delta\mathbf{u}_{i,j}d\Omega +\int_S\sigma_{ij}\mathbf{n}_j\delta\mathbf{u}_idS
$$

The stress tensor is symmetric $\sigma\_{ij}=\sigma\_{ji}$, hence:

$$
\int_\Omega\sigma_{ij,j}\delta\mathbf{u}_id\Omega=-\int_\Omega\sigma_{ij}\frac{1}{2}(\delta\mathbf{u}_{i,j}+\delta\mathbf{u}_{j,i})d\Omega +\int_S\sigma_{ij}\mathbf{n}_j\delta\mathbf{u}_idS
$$

where $ \frac{1}{2}(\delta\mathbf{u}_{i,j}+\delta\mathbf{u}_{j,i}) $ is the variational form of the Cauchy strain tensor $\epsilon = \frac{1}{2}(\mathbf{u}_{i,j}+\mathbf{u}_{j,i})$. Now the equation is formulated as:

$$
\int_\Omega\sigma_{ij}\delta\epsilon_{ij}d\Omega -\int_S\sigma_{ij}\mathbf{n}_j\delta\mathbf{u}_idS=0
$$

To be continued...

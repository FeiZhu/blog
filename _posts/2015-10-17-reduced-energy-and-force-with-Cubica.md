---
layout: post
title: Reduced Energy and Force with Cubica
---
This post presents some technical details about the representation of deformation energy and internal force in subspace using the method introduced in [Optimizing Cubature for Efficient Integration of Subspace Deformations](http://www.cs.cornell.edu/~djames/papers/cubature08.pdf).

Let's start with the basic notations:

* $\mathbf{X}$: the material coordinates in full space, $\mathbf{X}\in \mathbb{R}^{3\mathrm{n}}$.
* $\mathbf{x}$: the spatial coordinates in full space, $\mathbf{x}\in \mathbb{R}^{3\mathrm{n}}$.
* $\mathbf{q}$: the subspace deformation, $\mathbf{q}\in \mathbb{R}^{\mathrm{r}}$.
* $\mathbf{U}$: the displacement sensitivity matrix relating deformation in full space and subspace: $\mathbf{x} = \mathbf{X} + \mathbf{Uq}$, $\mathbf{U}\in \mathbb{R}^{3\mathrm{n}\times\mathrm{r}}$.

Now we take a step further to the deformation energy and internal force: 

* $\Psi(X,\mathbf{q})$ : the deformation energy density evaluated at material point $X\in \mathbb{R}^{3}$ with deformation represented in subspace $\mathbf{q}$ . For tetrahedral constant strain finite element,  $\Psi(X,\mathbf{q})$ is constant within the element.
* $W(\mathbf{q})$: the deformation energy evaluated with subspace deformation $\mathbf{q}$: $W(\mathbf{q}) = \int_\Omega\Psi(X,\mathbf{q})d\Omega_X$.
* $\mathbf{f}(\mathbf{q})$: the internal force evaluated with subspace deformation $\mathbf{q}$: 

$$
\begin{array}{l}
\qquad \mathbf{f}(\mathbf{q}) = -\int_\Omega\nabla_{\mathbf{q}}\Psi(X,\mathbf{q})d\Omega_X = \int_\mathbf{\Omega}\mathbf{g}(X,\mathbf{q})d\Omega_X
\end{array}
$$.

The reduced force density $\mathbf{g}(X,\mathbf{q})$ is the negative gradient of deformation energy density $\Psi(X,\mathbf{q})$ with respect to subspace deformation $ \mathbf{q}$ evaluated at material point $X$:

$$
\begin{array}{l}
\mathbf{g}(X,\mathbf{q}) &=-\nabla_\mathbf{q}\Psi(X,\mathbf{q})\\
                                        &= -\mathbf{U}^T\frac{\partial \Psi}{\partial \mathbf{x}}\\
\end{array}
$$

Entries of $\frac{\partial\Psi(X,\mathbf{q})}{\partial\mathbf{x}}$ vanish except those corresponding to vertices of the tetrahedral element that the material point $X$ lies in. Denote the element as $e$, from the [course note (P.29)](http://run.usc.edu/femdefo/barbic-courseNotes-modelReduction.pdf) we know that $[\frac{\partial \Psi}{\partial x_e^1}\  \frac{\partial \Psi}{\partial x_e^2}\  \frac{\partial \Psi}{\partial x_e^3}] = \mathbf{P}\mathbf{D}_m^{-T}$ and $  \frac{\partial \Psi}{\partial x_e^4} = -\frac{\partial \Psi}{\partial x_e^1}-\frac{\partial \Psi}{\partial x_e^2}-\frac{\partial \Psi}{\partial x_e^3}$.

In fact we do not need to explicitly assemble $\frac{\partial\Psi(X,\mathbf{q})}{\partial\mathbf{x}}$ to compute $\mathbf{g}(X,\mathbf{q}) $, the entries of $\mathbf{g}$ are computed as below:

$$
\mathbf{g}_\mathrm{i}(X,\mathbf{q}) = -\sum_{\mathbf{v}\in \mathbf{V}_e}\mathbf{U}_{\mathbf{v},\mathrm{i}}\frac{\partial\Psi}{\partial\mathbf{x}_\mathbf{v}},
$$

where $\mathbf{v} \in \mathbf{V} _{e}$
represents the vertices of the tetrahedral element $e$ that the material point $X$ lies in. To be more specific:

$$\mathbf{U}_{\mathbf{v},\mathrm{i}}\frac{\partial\Psi}{\partial\mathbf{x}_\mathbf{v}} = \sum_{a\in\{x,y,z\}}\mathbf{U}_{\mathbf{v}^a,\mathrm{i}}\frac{\partial\Psi}{\partial\mathbf{x}_{\mathbf{v}^a}}$$

and $a\in\{x,y,z\}$ is the three coordinate directions in Eulerian space.

As we've had the representation of reduced force density, it's trivial to express reduced internal force with Cubica approximation:

$$
\begin{array}{l}
\mathbf{f}(\mathbf{q}) & = \int_\mathbf{\Omega}\mathbf{g}(X,\mathbf{q})d\Omega_X\\
                                      &\approx \sum_{i=1}^N\omega_i\mathbf{g}(X_i,\mathbf{q}),
\end{array}
$$

where $N$ is the number of cubature elements and $\mathbf{g}(X_i,\mathbf{q})$ is evaluated at the $i$-th cubature element with reduced deformation $\mathbf{q}$.

It's easy to find that $\mathbf{g}(X,\mathbf{q})$ can be evaluated at $O(\mathbf{r})$ cost and therefore $\mathbf{f}(\mathbf{q})$ can be evaluated at $O(\mathbf{r}\times N)$ cost. If the number of cubature elements is at the scale of  $O(\mathbf{r})$, cubica method can evaluate internal forces at $O(\mathbf{r}^2)$ cost.

---
layout: post
title: Solve for Harmonic Displacements with FEM Discretization
---
For comparison reasons, I need to implement the SIGGRAPH2009 paper [Numerical Coarsening of Inhomogeneous Elastic Materials](http://users.cms.caltech.edu/~owhadi/publications/KMOD09.pdf). Everything else was set up except the harmonic displacements $ \mathbf{h}_{\alpha\beta} (1 \le \alpha \le \beta \le d) $:

$$
\begin{cases}
\textrm{div}(\mathbf{C}:\epsilon(\mathbf{h}_{\alpha\beta}))=0 \quad  &\textrm{inside}\  \Omega \\
(\mathbf{C}:\epsilon(\mathbf{h}_{\alpha\beta}))\cdot\mathbf{n}=\epsilon(x_{\alpha}\mathbf{e}_{\beta})\cdot\mathbf{n} \quad    &\textrm{for} \ \mathbf{x} \in \partial \Omega
\end{cases}
$$

This is in essence a static equilibrium problem with Neumann boundary conditions prescribing surface tractions equal to $ \epsilon(x\_{\alpha}\mathbf{e}\_{\beta})\cdot\mathbf{n} $. Here we consider solving this equation with FEM discretization.

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

where $\frac{1}{2}(\delta\mathbf{u}\_{i,j}+\delta\mathbf{u}\_{j,i})$ is the variational form of the Cauchy strain tensor $\epsilon = \frac{1}{2}(\mathbf{u}\_{i,j}+\mathbf{u}\_{j,i})$. The final variational form is:

$$
\int_\Omega\sigma_{ij}\delta\epsilon_{ij}d\Omega -\int_S\sigma_{ij}\mathbf{n}_j\delta\mathbf{u}_idS=0
$$

Before moving on, we introduce the [Voigt notation](http://en.wikipedia.org/wiki/Voigt_notation) where second-order tensors $ \sigma $ and $ \epsilon$ are expressed as 6x1 vectors, and fourth-order elastic tensor $ \mathbf{C} $ as a 6x6 matrix:

$$
\begin{array}{lcr}
[\sigma]  = [\sigma_{11},\sigma_{22},\sigma_{33},\sigma_{23},\sigma_{13},\sigma_{12}]^T \\
[\epsilon] =[\epsilon_{11},\epsilon_{22},\epsilon_{33},2\epsilon_{23},2\epsilon_{13},2\epsilon_{12}]^T\\
\end{array}
$$

With this convention, strain energy density $ \Psi = \sigma : \epsilon $ and the constitutive relation $\sigma = \mathbf{C}:\epsilon$ can be expressed as matrix productions $ \Psi = [\sigma]^T[\epsilon] $ and $ [\sigma] = [\mathbf{C}][\epsilon] $.

For Cauchy strains the relationship between strains and displacements is:

$$
[\epsilon]=[\mathbf{D}][\mathbf{u}]
$$

where $ [\mathbf{D}] $ is the matrix differentiation operator:

$$
[\mathbf{D}] =
\begin{bmatrix}
\frac{\partial}{\partial x} & 0 & 0 \\
0 & \frac{\partial}{\partial y} & 0 \\
0 & 0 & \frac{\partial}{\partial z} \\
0 & \frac{\partial}{\partial z} & \frac{\partial}{\partial y} \\
\frac{\partial}{\partial z} & 0 & \frac{\partial}{\partial x} \\
\frac{\partial}{\partial y} & \frac{\partial}{\partial x} & 0
\end{bmatrix}
$$

With FEM discretization, displacements at some point inside a finite element $ [\mathbf{u}] $ can be determined with the use of nodal displacements $ [\mathbf{q}] $ and shape functions $N\_i$ :

$$
\begin{array}{lcr}
[\mathbf{u}]=[\mathbf{N}][\mathbf{q}] \\
[\mathbf{N}] = 
\begin{bmatrix}
N_1 & 0 & 0 & N_2 & 0 & 0 & \dots & N_n & 0 & 0 \\
0 & N_1 & 0 & 0 & N_2 & 0 & \dots & 0 & N_n & 0 \\
0 & 0 & N_1 & 0 & 0 & N_2 & \dots & 0 & 0 & N_n
\end{bmatrix}
\end{array}
$$

Note: $ n $ is the number of vertices for each finite element.

Express the strains at some point inside a finite element with nodal displacements as:

$$
\begin{array}{lcr}
[\epsilon] = [\mathbf{B}][\mathbf{q}] \\
[\mathbf{B}] = [\mathbf{D}][\mathbf{N}] = [B_1 \ B_2 \ B_3 \dots \ B_n]
\end{array}
$$

where $ [\mathbf{B}] $ is the displacement differentiation matrix:

$$
[\mathbf{B}_i] =
\begin{bmatrix}
\frac{\partial N_i}{\partial x} & 0 & 0 \\
0 & \frac{\partial N_i}{\partial y} & 0 \\
0 & 0 & \frac{\partial N_i}{\partial z} \\
0 & \frac{\partial N_i}{\partial z} & \frac{\partial N_i}{\partial y} \\
\frac{\partial N_i}{\partial z} & 0 & \frac{\partial N_i}{\partial x} \\
\frac{\partial N_i}{\partial y} & \frac{\partial N_i}{\partial x} & 0
\end{bmatrix}
$$

Now let's go back to the variational  form. Using the strain-stress relations and the strain-displacements relations, we are able express the equation with nodal displacements:

$$
\begin{array}{lcr}
\quad\ \int_\Omega [\sigma]^T[\delta\epsilon]d\Omega-\int_S([\sigma][\mathbf{n}])^T[\delta\mathbf{u}]dS=0 \\
\Rightarrow \int_\Omega [\mathbf{q}]^T[\mathbf{B}]^T[\mathbf{C}]^T[\mathbf{B}][\delta\mathbf{q}]d\Omega-\int_S[\mathbf{t}]^T[\mathbf{N}][\delta\mathbf{q}]dS=0 \\
\Rightarrow \int_\Omega [\mathbf{q}]^T[\mathbf{B}]^T[\mathbf{C}]^T[\mathbf{B}]d\Omega-\int_S[\mathbf{t}]^T[\mathbf{N}]dS=0 \\
\Rightarrow \int_\Omega [\mathbf{B}]^T[\mathbf{C}][\mathbf{B}]d\Omega[\mathbf{q}]-\int_S[\mathbf{N}]^T[\mathbf{t}]dS=0 
\end{array}
$$

where $ \mathbf{t} $ is the traction at some point on surface. You may have noticed that this is the classical form of $ [\mathbf{k}][\mathbf{q}] = [\mathbf{f}]$:

$$
\begin{array}{clr}
[\mathbf{k}]=\int_\Omega [\mathbf{B}]^T[\mathbf{C}][\mathbf{B}]d\Omega \\
[\mathbf{f}]=\int_S[\mathbf{N}]^T[\mathbf{t}]dS
\end{array}
$$

Till now we've obtained the equation for each finite element $ [\mathbf{k}\_i][\mathbf{q}\_i] = [\mathbf{f}\_i] $, and we're going to assemble them to form the global equation system $ [\mathbf{K}][\mathbf{Q}] =[\mathbf{F}] $ . The process is usually called assembly.

Let us first consider directly assemble these element matrices in the simplest manner:

$$
\begin{array}{clr}
[\mathbf{Q}_d]=[[\mathbf{q}_1]^T \ \ [\mathbf{q}_2]^T \ \ \dots \ \ [\mathbf{q}_n]^T]^T \\
[\mathbf{F}_d]=[[\mathbf{f}_1]^T \ \ [\mathbf{f}_2]^T \ \ \dots \ \ [\mathbf{f}_n]^T]^T \\
[\mathbf{K}_d]=
\begin{bmatrix}
[\mathbf{k}_1] & 0 & 0 \\
0 & [\mathbf{k}_2] & 0 \\
0 & 0 & \dots
\end{bmatrix} 
\end{array}
$$

Note: here $ n $ is the number of finite elements.

It is evident that the directly assembled matrices are not what we want. But it is easy to find a matrix $[\mathbf{A}] $ such that:

$$
\begin{array}{clr}
[\mathbf{Q}_d]=[\mathbf{A}][\mathbf{Q}]\\
[\mathbf{F}_d]=[\mathbf{A}][\mathbf{F}]
\end{array}
$$

The matrix $ [\mathbf{A}] $ provides the transformation from global to local enumeration. $[\mathbf{A}\_{ij}] = 1$ if $ [\mathbf{Q}\_d]\_i $ corresponds to global vertex $ [\mathbf{Q}]\_j $, otherwise $[\mathbf{A}\_{ij}] = 0$ .

Conceptually speaking, the assembly process sums  up the contributions at a global vertex from the finite elements that share this common vertex.

In our case, the surface traction $ [\mathbf{f}] $ vanishes for internal element surfaces because  contributions from the two elements sharing that face are equal. While for external surfaces, the traction is prescribed boundary condition $ \epsilon(x\_{\alpha}\mathbf{e}\_{\beta})\cdot\mathbf{n} $.

To summarize, we briefly describe the procedures for implementation. First form the global stiffness matrix $[\mathbf{K}] $ using the local element stiffness matrices $[\mathbf{k}] $. For the right hand side, set the entry to zero if the vertex is an internal vertex. If the vertex belongs to some external faces, then for each external face, compute the traction at this vertex following the procedures below and sum up the contributions from all the faces:

1. compute the strain as $\epsilon(x\_\alpha\mathbf{e}\_\beta) = \frac{1}{2}(\mathbf{e}\_\alpha\otimes\mathbf{e}\_\beta+\mathbf{e}\_\beta\otimes\mathbf{e}\_\alpha)$
2.  compute the face normal $\mathbf{n}$ by cross producting the edges
3.  compute the traction $ \mathbf{t}=\epsilon(x\_\alpha\mathbf{e}\_\beta)\cdot\mathbf{n} $
4.  compute the forces at vertices as $ [\mathbf{f}]=\int\_S[\mathbf{N}]^T[\mathbf{t}]dS $

After solving $ [\mathbf{K}][\mathbf{Q}] =[\mathbf{F}] $ with boudnary condition $\mathbf{t}\_{\alpha\beta}$, we get harmonic displacement $\mathbf{h}\_{\alpha\beta}$. 6 boundary conditions lead to a set of 6 displacements.

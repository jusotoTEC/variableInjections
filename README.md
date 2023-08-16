# MATLAB Code from Numerical Experiments in Paper *"Optimal modeling of nonlinear systems: method of variable injections"*

## Author

* Anatoli Torokhti (https://people.unisa.edu.au/anatoli.torokhti) - Email: anatoli.torokhti@unisa.edu.au
* Pablo Soto-Quiros (https://www.tec.ac.cr/juan-pablo-soto-quiros) - Email: jusoto@tec.ac.cr

Anatoli Torokhti is an Associate Professor from the *University of South Australia* (https://www.unisa.edu.au/) in Mawson Lakes, SA, Australia

Pablo Soto-Quiros is an Associate Professor from the *Instituto Tecnológico de Costa Rica* (https://www.tec.ac.cr/) in Cartago, Costa Rica


## Description

* This repository contains the MATLAB code for numerical experiments presented in the paper "*Optimal modeling of nonlinear systems: method of variable injections*". 
* This paper has been submitted for publication in a scientific journal. 
* Our work addresses the development and justification of the new approach to the modeling of nonlinear systems. Let $\mathcal{F}$ be  an unknown  input-output map of the system with random input and output ${\bf y}$ and ${\bf x}$, respectively. It is assumed that  ${\bf y}$ and ${\bf x}$  are available and covariance matrices formed from ${\bf y}$ and ${\bf x}$ are known. We  determine a model of $\mathcal{F}$  so that an associated error is minimized. To this end,  the model $\mathcal{T}_p$ is constructed as a sum of $p+1$ particular parts, in the form $\mathcal{T}_p({\bf y})=\Sigma\;G_jH_jQ_j({\bf v}_j)$  where $G_j$ and $ H_j$, for $j=0,\ldots, p$,  are matrices to be determined, and ${\bf v}_j$, for $j=1,\ldots,p$, is a special  random vector  called the injection. We denote ${\bf v}_0={\bf y}$. Further, $Q_j$ is a special transform aimed to facilitate the numerical realization of  model $\mathcal{T}_p$. It is determined in a way allowing us to optimally determine $G_j$ and $H_j$ as a solution of $p+1$ separate error minimization problems which are simpler than the original minimization problem. The empirical determination of injections ${\bf v}_1,\ldots, {\bf v}_p$ is considered. The proposed method has several degrees of freedom to diminish the associated error. They are `degree' $p$ of $\mathcal{T}_p$, choice of matrices $G_0, H_0, \ldots,$ $G_p, H_p$, dimensions of   matrices $G_0, H_0, \ldots,$ $G_p, H_p$ and injections ${\bf v}_1,\ldots, {\bf v}_p$, respectively.  In particular, it is shown that a variation of the injections in their dimensionality and special forms allows us to increase  the accuracy  of  the proposed model $\ttt_p$. The proposed approach differs from known techniques by its ingredients mentioned above. Four numerical examples are provided. In the end, the open problem is formulated.


<p align="center"><img width="1200" src="https://github.com/jusotoTEC/multifiltering_transform/blob/main/img/img1.png"></p>

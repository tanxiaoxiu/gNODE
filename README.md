# gNODE
gNODE: gLV model-informed neural ordinary differential equations for modeling microbial community dynamics

## Framework of eNODEconstr
![](https://github.com/tanxiaoxiu/gNODE/blob/master/Framework.png)


## (A) Training module 
Longitudinal microbial abundance data are used to train a NeuralODE whose dynamics are defined by the gLV equations, enabling joint estimation of species-specific growth rates and interaction coefficients. 

## (B) Prediction module
The trained model supports ecological network inference, forecasting of community trajectories, and simulation of system responses under perturbation scenarios. The inferred ecological parameters are further used to compute driver scores, which quantify each species’ contribution to community structure and dynamics.
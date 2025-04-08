# gNODE
gNODE: gLV model-informed neural ordinary differential equations for modeling microbial community dynamics

## Framework of eNODEconstr
![](https://github.com/tanxiaoxiu/gNODE/blob/master/Framework.png)


## (A) Training module 
The input data consists of time-series abundance data for microbes. gNODE estimates the parameters of the gLV ecology-guided Neural Ordinary Differential Equations (NeuralODEs) from the input data. 

## (B) Prediction module
Based on the trained model, it can predict the interactions between microbes, the dynamic trajectories of microbes, and the trajectories after perturbations. Additionally, the estimated parameters can be used to calculate the driver score metric. The driver score quantifies the contribution of each species to the community.

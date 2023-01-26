# RobustGAMLSS

In this Github repository you can find implementations of different robust gamboostLSS Families-objects, which are more robust than the original ones for fitting GAMLSS with gradient boosting.

Using the R add-on package "mboost" and "gamboostLSS", you can specify different distributions via "families = robust_NEWFAMILY(rob=c)" inside of the "gamboostLSS" and "glmboostLSS" function. For the robust families itself, there is the additional argument "rob", where the user can specify the robustness constant "c". Higher values of "c" are less robust than smaller values of "c". We recommend to chose the robustness constant in depenence of the data and modelled distribution (this is explained in more detail in the article "Robust Gradient Boosting for Generalized Additive Models for Location, Scale and Shape" from Speller et al.).

Example for GAMLSS regression via gradient boosting:
model <- gamboostLSS(y ~ x1 + x2, data= data.frame(y,x1,x2), families = robust_GaussianLSS(rob=c))

E.g.:
robust_GaussianLSS(rob=c) , robust_GammaLSS(rob=c), GaussianLSS(), GammaLSS()


```r
# R code for running an example on a bodyfat data set with different covariates:

source{Robust_gamboostLSS_Families.R} # loading robust gamboostLSS families-objects for robust Gaussian and robust Gamma 
# file contains also required R-Package "mboost","gamboostLSS", and also only for Gamma distribution relevant packages: "gamlss.dist","EnvStats"
# file contains also additional functions within the families-objects, which might be relevant for further/own robust families

data("bodyfat", package = "TH.data") # load bodyfat data
# outcome: DEXfat
# all available covariates (age + waistcirc + hipcirc + elbowbreadth + kneebreadth + anthro3a + anthro3b + anthro3c + anthro4) are modelled by seperate linear base-learner, when using glmboostLSS (and not gamboostLSS)

c_value_Gaussian <- c_generate_Gaussian(bodyfat$DEXfat) # default for tau=0.05

c_value_Gamma <- c_generate_Gamma(bodyfat$DEXfat) # default for tau=0.05


stopping=400

set.seed(321)
glmLSS_robustGaussian <- glmboostLSS(DEXfat~. , method = "noncyclic" , control = boost_control(mstop=stopping) , families = robust_GaussianLSS(stabilization = "MAD",rob=c_value_Gaussian), data = bodyfat)
cvr_robustGaussian <- cvrisk(glmLSS_robustGaussian) # default method is 25-fold bootstrap
coef(glmLSS_robustGaussian[mstop(cvr_robustGaussian)] , off2int=TRUE , which="") #  coefficients of glmLSS_Gaussian at optimal stopping iteration for cvrisk

set.seed(321)
glmLSS_robustGamma <- glmboostLSS(DEXfat~. , method = "noncyclic" , control = boost_control(mstop=stopping) , families = robust_GaussianLSS(stabilization = "MAD",rob=c_value_Gamma), data = bodyfat)
cvr_robustGamma <- cvrisk(glmLSS_robustGamma) # default method is 25-fold bootstrap
coef(glmLSS_robustGamma[mstop(cvr_robustGamma)] , off2int=TRUE , which="")  # coefficients of glmLSS_Gamma at optimal stopping iteration for cvrisk

plot(glmLSS_robustGaussian)
plot(glmLSS_robustGamma)


# for specific covariates:
set.seed(321)
glmLSS_robustGaussian_2 <- glmboostLSS(DEXfat ~ age + waistcirc , method = "noncyclic" , control = boost_control(mstop=stopping) , families = robust_GaussianLSS(stabilization = "MAD",rob=c_value_Gaussian), data = bodyfat)
#cvr_robustGaussian_2 <- cvrisk(glmLSS_robustGaussian_2) # default method is 25-fold bootstrap
#coef(glmLSS_robustGaussian_2[mstop(cvr_robustGaussian_2)],off2int=TRUE)

set.seed(321)
glmLSS_robustGamma_2 <- glmboostLSS(DEXfat ~ age + waistcirc , method = "noncyclic" , control = boost_control(mstop=stopping) , families = robust_GaussianLSS(stabilization = "MAD",rob=c_value_Gamma), data = bodyfat)
#cvr_robustGamma_2 <- cvrisk(glmLSS_robustGamma_2) # default method is 25-fold bootstrap
#coef(glmLSS_robustGamma_2[mstop(cvr_robustGamma_2)],off2int=TRUE)

```

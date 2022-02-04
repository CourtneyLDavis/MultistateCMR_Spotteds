###############################################
###############################################
# Code for model implemented in:
# Davis, C.L., D.J. Muñoz, S.M. Amburgey, C.R. Dinsmore, E.W. Teitsworth, and D.A.W. Miller
# "Multistate mark-recapture model to estimate sex-specific dispersal rates and distances for a wetland breeding amphibian metapopulation"


# THIS HEAVILY BORROWS FROM CODE PUBLISHED IN:
# Worthington, H., R. McCrea, R. King and R. Griffiths. (2019).
# "Estimating abundance from multiple sampling capture-recapture data via a multi-state multi-period stopover model"
# The Annals of Applied Statistics, 13(4): 2043-2064.
# https://doi.org/10.1214/19-AOAS1264

###############################################
###############################################

# -------------------------------------------------------------------------------
# This script defines the function to unpack the parameter vector for the female
# model that estimates inter-annual movement as a function of Euclidean 
# distance between sites as well as temporary emigration.
# -------------------------------------------------------------------------------

# -------------------------------------------------------------------------------
# Inputs: param - vector of parameters
#               - [1] - log n.missed
#               - [2:6] - recruitment probabilities r[1]:r[T-1] on logistic scale
#               - [7] - survival probability (constant) on logit scale
#               - [8] - intercept for logistic regression on arrivals (constant)
#               - [9:14] - log gradient for logistic regression on arrivals
#               - [15:25] - initial discrete state probabilities (state 1) on logit scale
#               - [26:31] - capture probabilities on logit scale
#               - [32] - intercept for logistic regression on retention (constant) 
#               - [33:38] - log gradient for logistic regression on retention
#               - [39] - sigma.primary - scale parameter for rate of decay of between-period dispersal
#               - [40] - between-period site fidelity (probability of no movement, constant)
#               - [41] - omega1 - temporary emigration (probability of unobserved -> unobserved)
#               - [42] - omega2 - temporary emigration (probability of observered -> unobserved)
#               
#         n - number of observed individuals 
#         T - number of years
#         K - number of occasions in each year (vector length T)
# Outputs: list containing:
#          n.missed - number of missed individuals
#          N - total population size
#          r - recruitment probabilities
#          rstar - conditional recruitment probabilities
#          s - survival probabilities
#          betaintercept - intercept for logistic regression on arrivals
#          betagradient - gradient for logistic regression on arrivals
#          betalogit - logistic regression on arrivals
#          betastar - conditional arrival probabilities
#          omega1 -
#          omega2 -
#          alpha - initial discrete state probabilities
#          p - capture probabilities
#          sigma.primary - scale parameter that determines the rate of decay between-period dispersal probability as a function of distance
#          fid - movement probabilities between periods
#          phiintercept - intercept for logistic regression on retention
#          phigradient - gradient for logistic regression on retention
#          philogit - logistic regression on retention
#          phistar - conditional retention probabilities
#          gamma - transition matrices between periods
#          gammat - transition matrices within periods
#          pmat - observation matrices
#          pione - state probabilities between periods
#          pionet - state probabilities within periods
# -------------------------------------------------------------------------------


param.unpack <- function(param, n, T, K, distSq)  {
  
  # number of missed animals
  n.missed <- exp(param[1])
  
  # total population size
  N <- n+n.missed
  
  # recruitment probabilities
  glogitr <- param[2:6]
  r <- exp(glogitr)/(1+sum(exp(glogitr)))
  r[T] <- 1-sum(r[1:(T-1)]) 
  rstar <- rep(0, T)
  for (t in 1:T)  {
    rstar[t] <- r[t]/sum(r[t:T])
  }
  
  # survival probabilities
  logits <- rep(param[7],T-1)
  s <- 1/(1+exp(-logits))
  
  # arrival probabilities
  betaintercept <- rep(param[8],T)
  betagradient <- exp(param[9:14])
  betalogit <- NULL
  for (t in 1:T)  {
    betalogitt <- rep(0, K[t])
    for (k in 1:K[t])  {
      betalogitt[k] <- 1/(1+exp(-(betaintercept[t]+betagradient[t]*k)))
    }
    # must reach probability of 1 by final occasion within each year
    betalogitt <- betalogitt/betalogitt[K[t]]
    betalogit[[t]] <- betalogitt
  }
  
  # conditional arrival probabilities
  betastar <- NULL
  for (t in 1:T)  {
    betastar[[t]] <- rep(0, K[t])
    betastar[[t]][1] <- betalogit[[t]][1]
    for (k in 2:K[t])  {
      diff <- betalogit[[t]][k] - betalogit[[t]][k-1]
      if (isTRUE(diff > 0 | betalogit[[t]][k] < 0.5))  {
        betastar[[t]][k] <- (betalogit[[t]][k]-betalogit[[t]][k-1])/(1-betalogit[[t]][k-1])
      }
    }
  }

  
  # initial discrete state probabilities
  alpha <- matrix(0,nrow=12,ncol=1)
  alpha[1] <- exp(param[15])/(1 + exp(param[15]) + exp(param[16]) + exp(param[17]) + exp(param[18]) + exp(param[19]) + exp(param[20]) + exp(param[21]) + exp(param[22]) + exp(param[23]) + exp(param[24]) + exp(param[25]))
  alpha[2] <- exp(param[16])/(1 + exp(param[15]) + exp(param[16]) + exp(param[17]) + exp(param[18]) + exp(param[19]) + exp(param[20]) + exp(param[21]) + exp(param[22]) + exp(param[23]) + exp(param[24]) + exp(param[25]))
  alpha[3] <- exp(param[17])/(1 + exp(param[15]) + exp(param[16]) + exp(param[17]) + exp(param[18]) + exp(param[19]) + exp(param[20]) + exp(param[21]) + exp(param[22]) + exp(param[23]) + exp(param[24]) + exp(param[25]))
  alpha[4] <- exp(param[18])/(1 + exp(param[15]) + exp(param[16]) + exp(param[17]) + exp(param[18]) + exp(param[19]) + exp(param[20]) + exp(param[21]) + exp(param[22]) + exp(param[23]) + exp(param[24]) + exp(param[25]))
  alpha[5] <- exp(param[19])/(1 + exp(param[15]) + exp(param[16]) + exp(param[17]) + exp(param[18]) + exp(param[19]) + exp(param[20]) + exp(param[21]) + exp(param[22]) + exp(param[23]) + exp(param[24]) + exp(param[25]))
  alpha[6] <- exp(param[20])/(1 + exp(param[15]) + exp(param[16]) + exp(param[17]) + exp(param[18]) + exp(param[19]) + exp(param[20]) + exp(param[21]) + exp(param[22]) + exp(param[23]) + exp(param[24]) + exp(param[25]))
  alpha[7] <- exp(param[21])/(1 + exp(param[15]) + exp(param[16]) + exp(param[17]) + exp(param[18]) + exp(param[19]) + exp(param[20]) + exp(param[21]) + exp(param[22]) + exp(param[23]) + exp(param[24]) + exp(param[25]))
  alpha[8] <- exp(param[22])/(1 + exp(param[15]) + exp(param[16]) + exp(param[17]) + exp(param[18]) + exp(param[19]) + exp(param[20]) + exp(param[21]) + exp(param[22]) + exp(param[23]) + exp(param[24]) + exp(param[25]))
  alpha[9] <- exp(param[23])/(1 + exp(param[15]) + exp(param[16]) + exp(param[17]) + exp(param[18]) + exp(param[19]) + exp(param[20]) + exp(param[21]) + exp(param[22]) + exp(param[23]) + exp(param[24]) + exp(param[25]))
  alpha[10] <- exp(param[24])/(1 + exp(param[15]) + exp(param[16]) + exp(param[17]) + exp(param[18]) + exp(param[19]) + exp(param[20]) + exp(param[21]) + exp(param[22]) + exp(param[23]) + exp(param[24]) + exp(param[25]))
  alpha[11] <- exp(param[25])/(1 + exp(param[15]) + exp(param[16]) + exp(param[17]) + exp(param[18]) + exp(param[19]) + exp(param[20]) + exp(param[21]) + exp(param[22]) + exp(param[23]) + exp(param[24]) + exp(param[25]))
  alpha[12] <- 1-(alpha[1]+alpha[2]+alpha[3]+alpha[4]+alpha[5]+alpha[6]+alpha[7]+alpha[8]+alpha[9]+alpha[10]+alpha[11])
  
  
  # capture probabilities
  logitp <- param[26:31]
  pvec <- 1/(1+exp(-logitp))
  p <- vector('list', T)
  for (t in 1:T)  {
    p[[t]] <- pvec[1]
    pvec <- pvec[-1]
  }
  
  # retention probabilities
  phiintercept <- rep(param[32],T)
  phigradient <- -exp(param[33:38])
  philogit <- NULL
  for (t in 1:T)  {
    philogitt <- rep(0, K[t]-1)
    for (k in 1:(K[t]-1))  {
      philogitt[k] <- 1/(1+exp(-(phiintercept[t]+phigradient[t]*k)))
    }
    # must start at 1 on occasion 1
    philogitt <- philogitt/philogitt[1]
    if (!is.nan(philogitt[1]))  {
      philogit[[t]] <- philogitt
    } else if (is.nan(philogitt[1]))  {
      philogitt <- rep(0, K[t]-1)
      philogitt[1] <- 1
      philogit[[t]] <- philogitt
    }
  }
  
  # conditional retention probabilities
  phistar <- NULL
  for (t in 1:T)  {
    phistar[[t]] <- rep(0, K[t]-1)
    phistar[[t]][1] <- philogit[[t]][1]
    for (k in 2:(K[t]-1))  {
      diff <- philogit[[t]][k-1] - philogit[[t]][k]
      if (isTRUE(diff > 0 | philogit[[t]][k-1] > 0.5))  {
        phistar[[t]][k] <- philogit[[t]][k]/philogit[[t]][k-1]
      }
    }
  }
  
  # site fidelity matrix and movement between years  
  sigma.primary <- exp(param[39])
  Fid <- 1/(1+exp(-param[40]))
  
  # movement as a function of Euclidean distance between sites
  move.primary <- array(0, dim = c(12,12))
  move.primary[1,2] <- exp((-distSq[1,2]^2)/(2*sigma.primary^2))
  move.primary[1,3] <- exp((-distSq[1,3]^2)/(2*sigma.primary^2))
  move.primary[1,4] <- exp((-distSq[1,4]^2)/(2*sigma.primary^2))
  move.primary[1,5] <- exp((-distSq[1,5]^2)/(2*sigma.primary^2))
  move.primary[1,6] <- exp((-distSq[1,6]^2)/(2*sigma.primary^2))
  move.primary[1,7] <- exp((-distSq[1,7]^2)/(2*sigma.primary^2))
  move.primary[1,8] <- exp((-distSq[1,8]^2)/(2*sigma.primary^2))
  move.primary[1,9] <- exp((-distSq[1,9]^2)/(2*sigma.primary^2))
  move.primary[1,10] <- exp((-distSq[1,10]^2)/(2*sigma.primary^2))
  move.primary[1,11] <- exp((-distSq[1,11]^2)/(2*sigma.primary^2))
  move.primary[1,12] <- exp((-distSq[1,12]^2)/(2*sigma.primary^2))
  
  move.primary[2,1] <- exp((-distSq[2,1]^2)/(2*sigma.primary^2))
  move.primary[2,3] <- exp((-distSq[2,3]^2)/(2*sigma.primary^2))
  move.primary[2,4] <- exp((-distSq[2,4]^2)/(2*sigma.primary^2))
  move.primary[2,5] <- exp((-distSq[2,5]^2)/(2*sigma.primary^2))
  move.primary[2,6] <- exp((-distSq[2,6]^2)/(2*sigma.primary^2))
  move.primary[2,7] <- exp((-distSq[2,7]^2)/(2*sigma.primary^2))
  move.primary[2,8] <- exp((-distSq[2,8]^2)/(2*sigma.primary^2))
  move.primary[2,9] <- exp((-distSq[2,9]^2)/(2*sigma.primary^2))
  move.primary[2,10] <- exp((-distSq[2,10]^2)/(2*sigma.primary^2))
  move.primary[2,11] <- exp((-distSq[2,11]^2)/(2*sigma.primary^2))
  move.primary[2,12] <- exp((-distSq[2,12]^2)/(2*sigma.primary^2))
  
  move.primary[3,1] <- exp((-distSq[3,1]^2)/(2*sigma.primary^2))
  move.primary[3,2] <- exp((-distSq[3,2]^2)/(2*sigma.primary^2))
  move.primary[3,4] <- exp((-distSq[3,4]^2)/(2*sigma.primary^2))
  move.primary[3,5] <- exp((-distSq[3,5]^2)/(2*sigma.primary^2))
  move.primary[3,6] <- exp((-distSq[3,6]^2)/(2*sigma.primary^2))
  move.primary[3,7] <- exp((-distSq[3,7]^2)/(2*sigma.primary^2))
  move.primary[3,8] <- exp((-distSq[3,8]^2)/(2*sigma.primary^2))
  move.primary[3,9] <- exp((-distSq[3,9]^2)/(2*sigma.primary^2))
  move.primary[3,10] <- exp((-distSq[3,10]^2)/(2*sigma.primary^2))
  move.primary[3,11] <- exp((-distSq[3,11]^2)/(2*sigma.primary^2))
  move.primary[3,12] <- exp((-distSq[3,12]^2)/(2*sigma.primary^2))
  
  move.primary[4,1] <- exp((-distSq[4,1]^2)/(2*sigma.primary^2))
  move.primary[4,2] <- exp((-distSq[4,2]^2)/(2*sigma.primary^2))
  move.primary[4,3] <- exp((-distSq[4,3]^2)/(2*sigma.primary^2))
  move.primary[4,5] <- exp((-distSq[4,5]^2)/(2*sigma.primary^2))
  move.primary[4,6] <- exp((-distSq[4,6]^2)/(2*sigma.primary^2))
  move.primary[4,7] <- exp((-distSq[4,7]^2)/(2*sigma.primary^2))
  move.primary[4,8] <- exp((-distSq[4,8]^2)/(2*sigma.primary^2))
  move.primary[4,9] <- exp((-distSq[4,9]^2)/(2*sigma.primary^2))
  move.primary[4,10] <- exp((-distSq[4,10]^2)/(2*sigma.primary^2))
  move.primary[4,11] <- exp((-distSq[4,11]^2)/(2*sigma.primary^2))
  move.primary[4,12] <- exp((-distSq[4,12]^2)/(2*sigma.primary^2))
  
  move.primary[5,1] <- exp((-distSq[5,1]^2)/(2*sigma.primary^2))
  move.primary[5,2] <- exp((-distSq[5,2]^2)/(2*sigma.primary^2))
  move.primary[5,3] <- exp((-distSq[5,3]^2)/(2*sigma.primary^2))
  move.primary[5,4] <- exp((-distSq[5,4]^2)/(2*sigma.primary^2))
  move.primary[5,6] <- exp((-distSq[5,6]^2)/(2*sigma.primary^2))
  move.primary[5,7] <- exp((-distSq[5,7]^2)/(2*sigma.primary^2))
  move.primary[5,8] <- exp((-distSq[5,8]^2)/(2*sigma.primary^2))
  move.primary[5,9] <- exp((-distSq[5,9]^2)/(2*sigma.primary^2))
  move.primary[5,10] <- exp((-distSq[5,10]^2)/(2*sigma.primary^2))
  move.primary[5,11] <- exp((-distSq[5,11]^2)/(2*sigma.primary^2))
  move.primary[5,12] <- exp((-distSq[5,12]^2)/(2*sigma.primary^2))
  
  move.primary[6,1] <- exp((-distSq[6,1]^2)/(2*sigma.primary^2))
  move.primary[6,2] <- exp((-distSq[6,2]^2)/(2*sigma.primary^2))
  move.primary[6,3] <- exp((-distSq[6,3]^2)/(2*sigma.primary^2))
  move.primary[6,4] <- exp((-distSq[6,4]^2)/(2*sigma.primary^2))
  move.primary[6,5] <- exp((-distSq[6,5]^2)/(2*sigma.primary^2))
  move.primary[6,7] <- exp((-distSq[6,7]^2)/(2*sigma.primary^2))
  move.primary[6,8] <- exp((-distSq[6,8]^2)/(2*sigma.primary^2))
  move.primary[6,9] <- exp((-distSq[6,9]^2)/(2*sigma.primary^2))
  move.primary[6,10] <- exp((-distSq[6,10]^2)/(2*sigma.primary^2))
  move.primary[6,11] <- exp((-distSq[6,11]^2)/(2*sigma.primary^2))
  move.primary[6,12] <- exp((-distSq[6,12]^2)/(2*sigma.primary^2))
  
  move.primary[7,1] <- exp((-distSq[7,1]^2)/(2*sigma.primary^2))
  move.primary[7,2] <- exp((-distSq[7,2]^2)/(2*sigma.primary^2))
  move.primary[7,3] <- exp((-distSq[7,3]^2)/(2*sigma.primary^2))
  move.primary[7,4] <- exp((-distSq[7,4]^2)/(2*sigma.primary^2))
  move.primary[7,5] <- exp((-distSq[7,5]^2)/(2*sigma.primary^2))
  move.primary[7,6] <- exp((-distSq[7,6]^2)/(2*sigma.primary^2))
  move.primary[7,8] <- exp((-distSq[7,8]^2)/(2*sigma.primary^2))
  move.primary[7,9] <- exp((-distSq[7,9]^2)/(2*sigma.primary^2))
  move.primary[7,10] <- exp((-distSq[7,10]^2)/(2*sigma.primary^2))
  move.primary[7,11] <- exp((-distSq[7,11]^2)/(2*sigma.primary^2))
  move.primary[7,12] <- exp((-distSq[7,12]^2)/(2*sigma.primary^2))
  
  move.primary[8,1] <- exp((-distSq[8,1]^2)/(2*sigma.primary^2))
  move.primary[8,2] <- exp((-distSq[8,2]^2)/(2*sigma.primary^2))
  move.primary[8,3] <- exp((-distSq[8,3]^2)/(2*sigma.primary^2))
  move.primary[8,4] <- exp((-distSq[8,4]^2)/(2*sigma.primary^2))
  move.primary[8,5] <- exp((-distSq[8,5]^2)/(2*sigma.primary^2))
  move.primary[8,6] <- exp((-distSq[8,6]^2)/(2*sigma.primary^2))
  move.primary[8,7] <- exp((-distSq[8,7]^2)/(2*sigma.primary^2))
  move.primary[8,9] <- exp((-distSq[8,9]^2)/(2*sigma.primary^2))
  move.primary[8,10] <- exp((-distSq[8,10]^2)/(2*sigma.primary^2))
  move.primary[8,11] <- exp((-distSq[8,11]^2)/(2*sigma.primary^2))
  move.primary[8,12] <- exp((-distSq[8,12]^2)/(2*sigma.primary^2))
  
  move.primary[9,1] <- exp((-distSq[9,1]^2)/(2*sigma.primary^2))
  move.primary[9,2] <- exp((-distSq[9,2]^2)/(2*sigma.primary^2))
  move.primary[9,3] <- exp((-distSq[9,3]^2)/(2*sigma.primary^2))
  move.primary[9,4] <- exp((-distSq[9,4]^2)/(2*sigma.primary^2))
  move.primary[9,5] <- exp((-distSq[9,5]^2)/(2*sigma.primary^2))
  move.primary[9,6] <- exp((-distSq[9,6]^2)/(2*sigma.primary^2))
  move.primary[9,7] <- exp((-distSq[9,7]^2)/(2*sigma.primary^2))
  move.primary[9,8] <- exp((-distSq[9,8]^2)/(2*sigma.primary^2))
  move.primary[9,10] <- exp((-distSq[9,10]^2)/(2*sigma.primary^2))
  move.primary[9,11] <- exp((-distSq[9,11]^2)/(2*sigma.primary^2))
  move.primary[9,12] <- exp((-distSq[9,12]^2)/(2*sigma.primary^2))
  
  move.primary[10,1] <- exp((-distSq[10,1]^2)/(2*sigma.primary^2))
  move.primary[10,2] <- exp((-distSq[10,2]^2)/(2*sigma.primary^2))
  move.primary[10,3] <- exp((-distSq[10,3]^2)/(2*sigma.primary^2))
  move.primary[10,4] <- exp((-distSq[10,4]^2)/(2*sigma.primary^2))
  move.primary[10,5] <- exp((-distSq[10,5]^2)/(2*sigma.primary^2))
  move.primary[10,6] <- exp((-distSq[10,6]^2)/(2*sigma.primary^2))
  move.primary[10,7] <- exp((-distSq[10,7]^2)/(2*sigma.primary^2))
  move.primary[10,8] <- exp((-distSq[10,8]^2)/(2*sigma.primary^2))
  move.primary[10,9] <- exp((-distSq[10,9]^2)/(2*sigma.primary^2))
  move.primary[10,11] <- exp((-distSq[10,11]^2)/(2*sigma.primary^2))
  move.primary[10,12] <- exp((-distSq[10,12]^2)/(2*sigma.primary^2))
  
  move.primary[11,1] <- exp((-distSq[11,1]^2)/(2*sigma.primary^2))
  move.primary[11,2] <- exp((-distSq[11,2]^2)/(2*sigma.primary^2))
  move.primary[11,3] <- exp((-distSq[11,3]^2)/(2*sigma.primary^2))
  move.primary[11,4] <- exp((-distSq[11,4]^2)/(2*sigma.primary^2))
  move.primary[11,5] <- exp((-distSq[11,5]^2)/(2*sigma.primary^2))
  move.primary[11,6] <- exp((-distSq[11,6]^2)/(2*sigma.primary^2))
  move.primary[11,7] <- exp((-distSq[11,7]^2)/(2*sigma.primary^2))
  move.primary[11,8] <- exp((-distSq[11,8]^2)/(2*sigma.primary^2))
  move.primary[11,9] <- exp((-distSq[11,9]^2)/(2*sigma.primary^2))
  move.primary[11,10] <- exp((-distSq[11,10]^2)/(2*sigma.primary^2))
  move.primary[11,12] <- exp((-distSq[11,12]^2)/(2*sigma.primary^2))
  
  move.primary[12,1] <- exp((-distSq[12,1]^2)/(2*sigma.primary^2))
  move.primary[12,2] <- exp((-distSq[12,2]^2)/(2*sigma.primary^2))
  move.primary[12,3] <- exp((-distSq[12,3]^2)/(2*sigma.primary^2))
  move.primary[12,4] <- exp((-distSq[12,4]^2)/(2*sigma.primary^2))
  move.primary[12,5] <- exp((-distSq[12,5]^2)/(2*sigma.primary^2))
  move.primary[12,6] <- exp((-distSq[12,6]^2)/(2*sigma.primary^2))
  move.primary[12,7] <- exp((-distSq[12,7]^2)/(2*sigma.primary^2))
  move.primary[12,8] <- exp((-distSq[12,8]^2)/(2*sigma.primary^2))
  move.primary[12,9] <- exp((-distSq[12,9]^2)/(2*sigma.primary^2))
  move.primary[12,10] <- exp((-distSq[12,10]^2)/(2*sigma.primary^2))
  move.primary[12,11] <- exp((-distSq[12,11]^2)/(2*sigma.primary^2))
  
  fid <- array(0, dim = c(12,12))
  fid[1,2] <- (1-Fid)*(move.primary[1,2]/(move.primary[1,2]+move.primary[1,3]+move.primary[1,4]+move.primary[1,5]+move.primary[1,6]+move.primary[1,7]+move.primary[1,8]+move.primary[1,9]+move.primary[1,10]+move.primary[1,11]+move.primary[1,12]))
  fid[1,3] <- (1-Fid)*(move.primary[1,3]/(move.primary[1,2]+move.primary[1,3]+move.primary[1,4]+move.primary[1,5]+move.primary[1,6]+move.primary[1,7]+move.primary[1,8]+move.primary[1,9]+move.primary[1,10]+move.primary[1,11]+move.primary[1,12]))
  fid[1,4] <- (1-Fid)*(move.primary[1,4]/(move.primary[1,2]+move.primary[1,3]+move.primary[1,4]+move.primary[1,5]+move.primary[1,6]+move.primary[1,7]+move.primary[1,8]+move.primary[1,9]+move.primary[1,10]+move.primary[1,11]+move.primary[1,12]))
  fid[1,5] <- (1-Fid)*(move.primary[1,5]/(move.primary[1,2]+move.primary[1,3]+move.primary[1,4]+move.primary[1,5]+move.primary[1,6]+move.primary[1,7]+move.primary[1,8]+move.primary[1,9]+move.primary[1,10]+move.primary[1,11]+move.primary[1,12]))
  fid[1,6] <- (1-Fid)*(move.primary[1,6]/(move.primary[1,2]+move.primary[1,3]+move.primary[1,4]+move.primary[1,5]+move.primary[1,6]+move.primary[1,7]+move.primary[1,8]+move.primary[1,9]+move.primary[1,10]+move.primary[1,11]+move.primary[1,12]))
  fid[1,7] <- (1-Fid)*(move.primary[1,7]/(move.primary[1,2]+move.primary[1,3]+move.primary[1,4]+move.primary[1,5]+move.primary[1,6]+move.primary[1,7]+move.primary[1,8]+move.primary[1,9]+move.primary[1,10]+move.primary[1,11]+move.primary[1,12]))
  fid[1,8] <- (1-Fid)*(move.primary[1,8]/(move.primary[1,2]+move.primary[1,3]+move.primary[1,4]+move.primary[1,5]+move.primary[1,6]+move.primary[1,7]+move.primary[1,8]+move.primary[1,9]+move.primary[1,10]+move.primary[1,11]+move.primary[1,12]))
  fid[1,9] <- (1-Fid)*(move.primary[1,9]/(move.primary[1,2]+move.primary[1,3]+move.primary[1,4]+move.primary[1,5]+move.primary[1,6]+move.primary[1,7]+move.primary[1,8]+move.primary[1,9]+move.primary[1,10]+move.primary[1,11]+move.primary[1,12]))
  fid[1,10] <- (1-Fid)*(move.primary[1,10]/(move.primary[1,2]+move.primary[1,3]+move.primary[1,4]+move.primary[1,5]+move.primary[1,6]+move.primary[1,7]+move.primary[1,8]+move.primary[1,9]+move.primary[1,10]+move.primary[1,11]+move.primary[1,12]))
  fid[1,11] <- (1-Fid)*(move.primary[1,11]/(move.primary[1,2]+move.primary[1,3]+move.primary[1,4]+move.primary[1,5]+move.primary[1,6]+move.primary[1,7]+move.primary[1,8]+move.primary[1,9]+move.primary[1,10]+move.primary[1,11]+move.primary[1,12]))
  fid[1,12] <- (1-Fid)*(move.primary[1,12]/(move.primary[1,2]+move.primary[1,3]+move.primary[1,4]+move.primary[1,5]+move.primary[1,6]+move.primary[1,7]+move.primary[1,8]+move.primary[1,9]+move.primary[1,10]+move.primary[1,11]+move.primary[1,12]))
  fid[1,1] <- Fid
  
  fid[2,1] <- (1-Fid)*(move.primary[2,1]/(move.primary[2,1]+move.primary[2,3]+move.primary[2,4]+move.primary[2,5]+move.primary[2,6]+move.primary[2,7]+move.primary[2,8]+move.primary[2,9]+move.primary[2,10]+move.primary[2,11]+move.primary[2,12]))
  fid[2,3] <- (1-Fid)*(move.primary[2,3]/(move.primary[2,1]+move.primary[2,3]+move.primary[2,4]+move.primary[2,5]+move.primary[2,6]+move.primary[2,7]+move.primary[2,8]+move.primary[2,9]+move.primary[2,10]+move.primary[2,11]+move.primary[2,12]))
  fid[2,4] <- (1-Fid)*(move.primary[2,4]/(move.primary[2,1]+move.primary[2,3]+move.primary[2,4]+move.primary[2,5]+move.primary[2,6]+move.primary[2,7]+move.primary[2,8]+move.primary[2,9]+move.primary[2,10]+move.primary[2,11]+move.primary[2,12]))
  fid[2,5] <- (1-Fid)*(move.primary[2,5]/(move.primary[2,1]+move.primary[2,3]+move.primary[2,4]+move.primary[2,5]+move.primary[2,6]+move.primary[2,7]+move.primary[2,8]+move.primary[2,9]+move.primary[2,10]+move.primary[2,11]+move.primary[2,12]))
  fid[2,6] <- (1-Fid)*(move.primary[2,6]/(move.primary[2,1]+move.primary[2,3]+move.primary[2,4]+move.primary[2,5]+move.primary[2,6]+move.primary[2,7]+move.primary[2,8]+move.primary[2,9]+move.primary[2,10]+move.primary[2,11]+move.primary[2,12]))
  fid[2,7] <- (1-Fid)*(move.primary[2,7]/(move.primary[2,1]+move.primary[2,3]+move.primary[2,4]+move.primary[2,5]+move.primary[2,6]+move.primary[2,7]+move.primary[2,8]+move.primary[2,9]+move.primary[2,10]+move.primary[2,11]+move.primary[2,12]))
  fid[2,8] <- (1-Fid)*(move.primary[2,8]/(move.primary[2,1]+move.primary[2,3]+move.primary[2,4]+move.primary[2,5]+move.primary[2,6]+move.primary[2,7]+move.primary[2,8]+move.primary[2,9]+move.primary[2,10]+move.primary[2,11]+move.primary[2,12]))
  fid[2,9] <- (1-Fid)*(move.primary[2,9]/(move.primary[2,1]+move.primary[2,3]+move.primary[2,4]+move.primary[2,5]+move.primary[2,6]+move.primary[2,7]+move.primary[2,8]+move.primary[2,9]+move.primary[2,10]+move.primary[2,11]+move.primary[2,12]))
  fid[2,10] <- (1-Fid)*(move.primary[2,10]/(move.primary[2,1]+move.primary[2,3]+move.primary[2,4]+move.primary[2,5]+move.primary[2,6]+move.primary[2,7]+move.primary[2,8]+move.primary[2,9]+move.primary[2,10]+move.primary[2,11]+move.primary[2,12]))
  fid[2,11] <- (1-Fid)*(move.primary[2,11]/(move.primary[2,1]+move.primary[2,3]+move.primary[2,4]+move.primary[2,5]+move.primary[2,6]+move.primary[2,7]+move.primary[2,8]+move.primary[2,9]+move.primary[2,10]+move.primary[2,11]+move.primary[2,12]))
  fid[2,12] <- (1-Fid)*(move.primary[2,12]/(move.primary[2,1]+move.primary[2,3]+move.primary[2,4]+move.primary[2,5]+move.primary[2,6]+move.primary[2,7]+move.primary[2,8]+move.primary[2,9]+move.primary[2,10]+move.primary[2,11]+move.primary[2,12]))
  fid[2,2] <- Fid
  
  fid[3,1] <- (1-Fid)*(move.primary[3,1]/(move.primary[3,1]+move.primary[3,2]+move.primary[3,4]+move.primary[3,5]+move.primary[3,6]+move.primary[3,7]+move.primary[3,8]+move.primary[3,9]+move.primary[3,10]+move.primary[3,11]+move.primary[3,12]))
  fid[3,2] <- (1-Fid)*(move.primary[3,2]/(move.primary[3,1]+move.primary[3,2]+move.primary[3,4]+move.primary[3,5]+move.primary[3,6]+move.primary[3,7]+move.primary[3,8]+move.primary[3,9]+move.primary[3,10]+move.primary[3,11]+move.primary[3,12]))
  fid[3,4] <- (1-Fid)*(move.primary[3,4]/(move.primary[3,1]+move.primary[3,2]+move.primary[3,4]+move.primary[3,5]+move.primary[3,6]+move.primary[3,7]+move.primary[3,8]+move.primary[3,9]+move.primary[3,10]+move.primary[3,11]+move.primary[3,12]))
  fid[3,5] <- (1-Fid)*(move.primary[3,5]/(move.primary[3,1]+move.primary[3,2]+move.primary[3,4]+move.primary[3,5]+move.primary[3,6]+move.primary[3,7]+move.primary[3,8]+move.primary[3,9]+move.primary[3,10]+move.primary[3,11]+move.primary[3,12]))
  fid[3,6] <- (1-Fid)*(move.primary[3,6]/(move.primary[3,1]+move.primary[3,2]+move.primary[3,4]+move.primary[3,5]+move.primary[3,6]+move.primary[3,7]+move.primary[3,8]+move.primary[3,9]+move.primary[3,10]+move.primary[3,11]+move.primary[3,12]))
  fid[3,7] <- (1-Fid)*(move.primary[3,7]/(move.primary[3,1]+move.primary[3,2]+move.primary[3,4]+move.primary[3,5]+move.primary[3,6]+move.primary[3,7]+move.primary[3,8]+move.primary[3,9]+move.primary[3,10]+move.primary[3,11]+move.primary[3,12]))
  fid[3,8] <- (1-Fid)*(move.primary[3,8]/(move.primary[3,1]+move.primary[3,2]+move.primary[3,4]+move.primary[3,5]+move.primary[3,6]+move.primary[3,7]+move.primary[3,8]+move.primary[3,9]+move.primary[3,10]+move.primary[3,11]+move.primary[3,12]))
  fid[3,9] <- (1-Fid)*(move.primary[3,9]/(move.primary[3,1]+move.primary[3,2]+move.primary[3,4]+move.primary[3,5]+move.primary[3,6]+move.primary[3,7]+move.primary[3,8]+move.primary[3,9]+move.primary[3,10]+move.primary[3,11]+move.primary[3,12]))
  fid[3,10] <- (1-Fid)*(move.primary[3,10]/(move.primary[3,1]+move.primary[3,2]+move.primary[3,4]+move.primary[3,5]+move.primary[3,6]+move.primary[3,7]+move.primary[3,8]+move.primary[3,9]+move.primary[3,10]+move.primary[3,11]+move.primary[3,12]))
  fid[3,11] <- (1-Fid)*(move.primary[3,11]/(move.primary[3,1]+move.primary[3,2]+move.primary[3,4]+move.primary[3,5]+move.primary[3,6]+move.primary[3,7]+move.primary[3,8]+move.primary[3,9]+move.primary[3,10]+move.primary[3,11]+move.primary[3,12]))
  fid[3,12] <- (1-Fid)*(move.primary[3,12]/(move.primary[3,1]+move.primary[3,2]+move.primary[3,4]+move.primary[3,5]+move.primary[3,6]+move.primary[3,7]+move.primary[3,8]+move.primary[3,9]+move.primary[3,10]+move.primary[3,11]+move.primary[3,12]))
  fid[3,3] <- Fid
  
  fid[4,1] <- (1-Fid)*(move.primary[4,1]/(move.primary[4,1]+move.primary[4,2]+move.primary[4,3]+move.primary[4,5]+move.primary[4,6]+move.primary[4,7]+move.primary[4,8]+move.primary[4,9]+move.primary[4,10]+move.primary[4,11]+move.primary[4,12]))
  fid[4,2] <- (1-Fid)*(move.primary[4,2]/(move.primary[4,1]+move.primary[4,2]+move.primary[4,3]+move.primary[4,5]+move.primary[4,6]+move.primary[4,7]+move.primary[4,8]+move.primary[4,9]+move.primary[4,10]+move.primary[4,11]+move.primary[4,12]))
  fid[4,3] <- (1-Fid)*(move.primary[4,3]/(move.primary[4,1]+move.primary[4,2]+move.primary[4,3]+move.primary[4,5]+move.primary[4,6]+move.primary[4,7]+move.primary[4,8]+move.primary[4,9]+move.primary[4,10]+move.primary[4,11]+move.primary[4,12]))
  fid[4,5] <- (1-Fid)*(move.primary[4,5]/(move.primary[4,1]+move.primary[4,2]+move.primary[4,3]+move.primary[4,5]+move.primary[4,6]+move.primary[4,7]+move.primary[4,8]+move.primary[4,9]+move.primary[4,10]+move.primary[4,11]+move.primary[4,12]))
  fid[4,6] <- (1-Fid)*(move.primary[4,6]/(move.primary[4,1]+move.primary[4,2]+move.primary[4,3]+move.primary[4,5]+move.primary[4,6]+move.primary[4,7]+move.primary[4,8]+move.primary[4,9]+move.primary[4,10]+move.primary[4,11]+move.primary[4,12]))
  fid[4,7] <- (1-Fid)*(move.primary[4,7]/(move.primary[4,1]+move.primary[4,2]+move.primary[4,3]+move.primary[4,5]+move.primary[4,6]+move.primary[4,7]+move.primary[4,8]+move.primary[4,9]+move.primary[4,10]+move.primary[4,11]+move.primary[4,12]))
  fid[4,8] <- (1-Fid)*(move.primary[4,8]/(move.primary[4,1]+move.primary[4,2]+move.primary[4,3]+move.primary[4,5]+move.primary[4,6]+move.primary[4,7]+move.primary[4,8]+move.primary[4,9]+move.primary[4,10]+move.primary[4,11]+move.primary[4,12]))
  fid[4,9] <- (1-Fid)*(move.primary[4,9]/(move.primary[4,1]+move.primary[4,2]+move.primary[4,3]+move.primary[4,5]+move.primary[4,6]+move.primary[4,7]+move.primary[4,8]+move.primary[4,9]+move.primary[4,10]+move.primary[4,11]+move.primary[4,12]))
  fid[4,10] <- (1-Fid)*(move.primary[4,10]/(move.primary[4,1]+move.primary[4,2]+move.primary[4,3]+move.primary[4,5]+move.primary[4,6]+move.primary[4,7]+move.primary[4,8]+move.primary[4,9]+move.primary[4,10]+move.primary[4,11]+move.primary[4,12]))
  fid[4,11] <- (1-Fid)*(move.primary[4,11]/(move.primary[4,1]+move.primary[4,2]+move.primary[4,3]+move.primary[4,5]+move.primary[4,6]+move.primary[4,7]+move.primary[4,8]+move.primary[4,9]+move.primary[4,10]+move.primary[4,11]+move.primary[4,12]))
  fid[4,12] <- (1-Fid)*(move.primary[4,12]/(move.primary[4,1]+move.primary[4,2]+move.primary[4,3]+move.primary[4,5]+move.primary[4,6]+move.primary[4,7]+move.primary[4,8]+move.primary[4,9]+move.primary[4,10]+move.primary[4,11]+move.primary[4,12]))
  fid[4,4] <- Fid
  
  fid[5,1] <- (1-Fid)*(move.primary[5,1]/(move.primary[5,1]+move.primary[5,2]+move.primary[5,3]+move.primary[5,4]+move.primary[5,6]+move.primary[5,7]+move.primary[5,8]+move.primary[5,9]+move.primary[5,10]+move.primary[5,11]+move.primary[5,12]))
  fid[5,2] <- (1-Fid)*(move.primary[5,2]/(move.primary[5,1]+move.primary[5,2]+move.primary[5,3]+move.primary[5,4]+move.primary[5,6]+move.primary[5,7]+move.primary[5,8]+move.primary[5,9]+move.primary[5,10]+move.primary[5,11]+move.primary[5,12]))
  fid[5,3] <- (1-Fid)*(move.primary[5,3]/(move.primary[5,1]+move.primary[5,2]+move.primary[5,3]+move.primary[5,4]+move.primary[5,6]+move.primary[5,7]+move.primary[5,8]+move.primary[5,9]+move.primary[5,10]+move.primary[5,11]+move.primary[5,12]))
  fid[5,4] <- (1-Fid)*(move.primary[5,4]/(move.primary[5,1]+move.primary[5,2]+move.primary[5,3]+move.primary[5,4]+move.primary[5,6]+move.primary[5,7]+move.primary[5,8]+move.primary[5,9]+move.primary[5,10]+move.primary[5,11]+move.primary[5,12]))
  fid[5,6] <- (1-Fid)*(move.primary[5,6]/(move.primary[5,1]+move.primary[5,2]+move.primary[5,3]+move.primary[5,4]+move.primary[5,6]+move.primary[5,7]+move.primary[5,8]+move.primary[5,9]+move.primary[5,10]+move.primary[5,11]+move.primary[5,12]))
  fid[5,7] <- (1-Fid)*(move.primary[5,7]/(move.primary[5,1]+move.primary[5,2]+move.primary[5,3]+move.primary[5,4]+move.primary[5,6]+move.primary[5,7]+move.primary[5,8]+move.primary[5,9]+move.primary[5,10]+move.primary[5,11]+move.primary[5,12]))
  fid[5,8] <- (1-Fid)*(move.primary[5,8]/(move.primary[5,1]+move.primary[5,2]+move.primary[5,3]+move.primary[5,4]+move.primary[5,6]+move.primary[5,7]+move.primary[5,8]+move.primary[5,9]+move.primary[5,10]+move.primary[5,11]+move.primary[5,12]))
  fid[5,9] <- (1-Fid)*(move.primary[5,9]/(move.primary[5,1]+move.primary[5,2]+move.primary[5,3]+move.primary[5,4]+move.primary[5,6]+move.primary[5,7]+move.primary[5,8]+move.primary[5,9]+move.primary[5,10]+move.primary[5,11]+move.primary[5,12]))
  fid[5,10] <- (1-Fid)*(move.primary[5,10]/(move.primary[5,1]+move.primary[5,2]+move.primary[5,3]+move.primary[5,4]+move.primary[5,6]+move.primary[5,7]+move.primary[5,8]+move.primary[5,9]+move.primary[5,10]+move.primary[5,11]+move.primary[5,12]))
  fid[5,11] <- (1-Fid)*(move.primary[5,11]/(move.primary[5,1]+move.primary[5,2]+move.primary[5,3]+move.primary[5,4]+move.primary[5,6]+move.primary[5,7]+move.primary[5,8]+move.primary[5,9]+move.primary[5,10]+move.primary[5,11]+move.primary[5,12]))
  fid[5,12] <- (1-Fid)*(move.primary[5,12]/(move.primary[5,1]+move.primary[5,2]+move.primary[5,3]+move.primary[5,4]+move.primary[5,6]+move.primary[5,7]+move.primary[5,8]+move.primary[5,9]+move.primary[5,10]+move.primary[5,11]+move.primary[5,12]))
  fid[5,5] <- Fid
  
  fid[6,1] <- (1-Fid)*(move.primary[6,1]/(move.primary[6,1]+move.primary[6,2]+move.primary[6,3]+move.primary[6,4]+move.primary[6,5]+move.primary[6,7]+move.primary[6,8]+move.primary[6,9]+move.primary[6,10]+move.primary[6,11]+move.primary[6,12]))
  fid[6,2] <- (1-Fid)*(move.primary[6,2]/(move.primary[6,1]+move.primary[6,2]+move.primary[6,3]+move.primary[6,4]+move.primary[6,5]+move.primary[6,7]+move.primary[6,8]+move.primary[6,9]+move.primary[6,10]+move.primary[6,11]+move.primary[6,12]))
  fid[6,3] <- (1-Fid)*(move.primary[6,3]/(move.primary[6,1]+move.primary[6,2]+move.primary[6,3]+move.primary[6,4]+move.primary[6,5]+move.primary[6,7]+move.primary[6,8]+move.primary[6,9]+move.primary[6,10]+move.primary[6,11]+move.primary[6,12]))
  fid[6,4] <- (1-Fid)*(move.primary[6,4]/(move.primary[6,1]+move.primary[6,2]+move.primary[6,3]+move.primary[6,4]+move.primary[6,5]+move.primary[6,7]+move.primary[6,8]+move.primary[6,9]+move.primary[6,10]+move.primary[6,11]+move.primary[6,12]))
  fid[6,5] <- (1-Fid)*(move.primary[6,5]/(move.primary[6,1]+move.primary[6,2]+move.primary[6,3]+move.primary[6,4]+move.primary[6,5]+move.primary[6,7]+move.primary[6,8]+move.primary[6,9]+move.primary[6,10]+move.primary[6,11]+move.primary[6,12]))
  fid[6,7] <- (1-Fid)*(move.primary[6,7]/(move.primary[6,1]+move.primary[6,2]+move.primary[6,3]+move.primary[6,4]+move.primary[6,5]+move.primary[6,7]+move.primary[6,8]+move.primary[6,9]+move.primary[6,10]+move.primary[6,11]+move.primary[6,12]))
  fid[6,8] <- (1-Fid)*(move.primary[6,8]/(move.primary[6,1]+move.primary[6,2]+move.primary[6,3]+move.primary[6,4]+move.primary[6,5]+move.primary[6,7]+move.primary[6,8]+move.primary[6,9]+move.primary[6,10]+move.primary[6,11]+move.primary[6,12]))
  fid[6,9] <- (1-Fid)*(move.primary[6,9]/(move.primary[6,1]+move.primary[6,2]+move.primary[6,3]+move.primary[6,4]+move.primary[6,5]+move.primary[6,7]+move.primary[6,8]+move.primary[6,9]+move.primary[6,10]+move.primary[6,11]+move.primary[6,12]))
  fid[6,10] <- (1-Fid)*(move.primary[6,10]/(move.primary[6,1]+move.primary[6,2]+move.primary[6,3]+move.primary[6,4]+move.primary[6,5]+move.primary[6,7]+move.primary[6,8]+move.primary[6,9]+move.primary[6,10]+move.primary[6,11]+move.primary[6,12]))
  fid[6,11] <- (1-Fid)*(move.primary[6,11]/(move.primary[6,1]+move.primary[6,2]+move.primary[6,3]+move.primary[6,4]+move.primary[6,5]+move.primary[6,7]+move.primary[6,8]+move.primary[6,9]+move.primary[6,10]+move.primary[6,11]+move.primary[6,12]))
  fid[6,12] <- (1-Fid)*(move.primary[6,12]/(move.primary[6,1]+move.primary[6,2]+move.primary[6,3]+move.primary[6,4]+move.primary[6,5]+move.primary[6,7]+move.primary[6,8]+move.primary[6,9]+move.primary[6,10]+move.primary[6,11]+move.primary[6,12]))
  fid[6,6] <- Fid
  
  fid[7,1] <- (1-Fid)*(move.primary[7,1]/(move.primary[7,1]+move.primary[7,2]+move.primary[7,3]+move.primary[7,4]+move.primary[7,5]+move.primary[7,6]+move.primary[7,8]+move.primary[7,9]+move.primary[7,10]+move.primary[7,11]+move.primary[7,12]))
  fid[7,2] <- (1-Fid)*(move.primary[7,2]/(move.primary[7,1]+move.primary[7,2]+move.primary[7,3]+move.primary[7,4]+move.primary[7,5]+move.primary[7,6]+move.primary[7,8]+move.primary[7,9]+move.primary[7,10]+move.primary[7,11]+move.primary[7,12]))
  fid[7,3] <- (1-Fid)*(move.primary[7,3]/(move.primary[7,1]+move.primary[7,2]+move.primary[7,3]+move.primary[7,4]+move.primary[7,5]+move.primary[7,6]+move.primary[7,8]+move.primary[7,9]+move.primary[7,10]+move.primary[7,11]+move.primary[7,12]))
  fid[7,4] <- (1-Fid)*(move.primary[7,4]/(move.primary[7,1]+move.primary[7,2]+move.primary[7,3]+move.primary[7,4]+move.primary[7,5]+move.primary[7,6]+move.primary[7,8]+move.primary[7,9]+move.primary[7,10]+move.primary[7,11]+move.primary[7,12]))
  fid[7,5] <- (1-Fid)*(move.primary[7,5]/(move.primary[7,1]+move.primary[7,2]+move.primary[7,3]+move.primary[7,4]+move.primary[7,5]+move.primary[7,6]+move.primary[7,8]+move.primary[7,9]+move.primary[7,10]+move.primary[7,11]+move.primary[7,12]))
  fid[7,6] <- (1-Fid)*(move.primary[7,6]/(move.primary[7,1]+move.primary[7,2]+move.primary[7,3]+move.primary[7,4]+move.primary[7,5]+move.primary[7,6]+move.primary[7,8]+move.primary[7,9]+move.primary[7,10]+move.primary[7,11]+move.primary[7,12]))
  fid[7,8] <- (1-Fid)*(move.primary[7,8]/(move.primary[7,1]+move.primary[7,2]+move.primary[7,3]+move.primary[7,4]+move.primary[7,5]+move.primary[7,6]+move.primary[7,8]+move.primary[7,9]+move.primary[7,10]+move.primary[7,11]+move.primary[7,12]))
  fid[7,9] <- (1-Fid)*(move.primary[7,9]/(move.primary[7,1]+move.primary[7,2]+move.primary[7,3]+move.primary[7,4]+move.primary[7,5]+move.primary[7,6]+move.primary[7,8]+move.primary[7,9]+move.primary[7,10]+move.primary[7,11]+move.primary[7,12]))
  fid[7,10] <- (1-Fid)*(move.primary[7,10]/(move.primary[7,1]+move.primary[7,2]+move.primary[7,3]+move.primary[7,4]+move.primary[7,5]+move.primary[7,6]+move.primary[7,8]+move.primary[7,9]+move.primary[7,10]+move.primary[7,11]+move.primary[7,12]))
  fid[7,11] <- (1-Fid)*(move.primary[7,11]/(move.primary[7,1]+move.primary[7,2]+move.primary[7,3]+move.primary[7,4]+move.primary[7,5]+move.primary[7,6]+move.primary[7,8]+move.primary[7,9]+move.primary[7,10]+move.primary[7,11]+move.primary[7,12]))
  fid[7,12] <- (1-Fid)*(move.primary[7,12]/(move.primary[7,1]+move.primary[7,2]+move.primary[7,3]+move.primary[7,4]+move.primary[7,5]+move.primary[7,6]+move.primary[7,8]+move.primary[7,9]+move.primary[7,10]+move.primary[7,11]+move.primary[7,12]))
  fid[7,7] <- Fid
  
  fid[8,1] <- (1-Fid)*(move.primary[8,1]/(move.primary[8,1]+move.primary[8,2]+move.primary[8,3]+move.primary[8,4]+move.primary[8,5]+move.primary[8,6]+move.primary[8,7]+move.primary[8,9]+move.primary[8,10]+move.primary[8,11]+move.primary[8,12]))
  fid[8,2] <- (1-Fid)*(move.primary[8,2]/(move.primary[8,1]+move.primary[8,2]+move.primary[8,3]+move.primary[8,4]+move.primary[8,5]+move.primary[8,6]+move.primary[8,7]+move.primary[8,9]+move.primary[8,10]+move.primary[8,11]+move.primary[8,12]))
  fid[8,3] <- (1-Fid)*(move.primary[8,3]/(move.primary[8,1]+move.primary[8,2]+move.primary[8,3]+move.primary[8,4]+move.primary[8,5]+move.primary[8,6]+move.primary[8,7]+move.primary[8,9]+move.primary[8,10]+move.primary[8,11]+move.primary[8,12]))
  fid[8,4] <- (1-Fid)*(move.primary[8,4]/(move.primary[8,1]+move.primary[8,2]+move.primary[8,3]+move.primary[8,4]+move.primary[8,5]+move.primary[8,6]+move.primary[8,7]+move.primary[8,9]+move.primary[8,10]+move.primary[8,11]+move.primary[8,12]))
  fid[8,5] <- (1-Fid)*(move.primary[8,5]/(move.primary[8,1]+move.primary[8,2]+move.primary[8,3]+move.primary[8,4]+move.primary[8,5]+move.primary[8,6]+move.primary[8,7]+move.primary[8,9]+move.primary[8,10]+move.primary[8,11]+move.primary[8,12]))
  fid[8,6] <- (1-Fid)*(move.primary[8,6]/(move.primary[8,1]+move.primary[8,2]+move.primary[8,3]+move.primary[8,4]+move.primary[8,5]+move.primary[8,6]+move.primary[8,7]+move.primary[8,9]+move.primary[8,10]+move.primary[8,11]+move.primary[8,12]))
  fid[8,7] <- (1-Fid)*(move.primary[8,7]/(move.primary[8,1]+move.primary[8,2]+move.primary[8,3]+move.primary[8,4]+move.primary[8,5]+move.primary[8,6]+move.primary[8,7]+move.primary[8,9]+move.primary[8,10]+move.primary[8,11]+move.primary[8,12]))
  fid[8,9] <- (1-Fid)*(move.primary[8,9]/(move.primary[8,1]+move.primary[8,2]+move.primary[8,3]+move.primary[8,4]+move.primary[8,5]+move.primary[8,6]+move.primary[8,7]+move.primary[8,9]+move.primary[8,10]+move.primary[8,11]+move.primary[8,12]))
  fid[8,10] <- (1-Fid)*(move.primary[8,10]/(move.primary[8,1]+move.primary[8,2]+move.primary[8,3]+move.primary[8,4]+move.primary[8,5]+move.primary[8,6]+move.primary[8,7]+move.primary[8,9]+move.primary[8,10]+move.primary[8,11]+move.primary[8,12]))
  fid[8,11] <- (1-Fid)*(move.primary[8,11]/(move.primary[8,1]+move.primary[8,2]+move.primary[8,3]+move.primary[8,4]+move.primary[8,5]+move.primary[8,6]+move.primary[8,7]+move.primary[8,9]+move.primary[8,10]+move.primary[8,11]+move.primary[8,12]))
  fid[8,12] <- (1-Fid)*(move.primary[8,12]/(move.primary[8,1]+move.primary[8,2]+move.primary[8,3]+move.primary[8,4]+move.primary[8,5]+move.primary[8,6]+move.primary[8,7]+move.primary[8,9]+move.primary[8,10]+move.primary[8,11]+move.primary[8,12]))
  fid[8,8] <- Fid
  
  fid[9,1] <- (1-Fid)*(move.primary[9,1]/(move.primary[9,1]+move.primary[9,2]+move.primary[9,3]+move.primary[9,4]+move.primary[9,5]+move.primary[9,6]+move.primary[9,7]+move.primary[9,8]+move.primary[9,10]+move.primary[9,11]+move.primary[9,12]))
  fid[9,2] <- (1-Fid)*(move.primary[9,2]/(move.primary[9,1]+move.primary[9,2]+move.primary[9,3]+move.primary[9,4]+move.primary[9,5]+move.primary[9,6]+move.primary[9,7]+move.primary[9,8]+move.primary[9,10]+move.primary[9,11]+move.primary[9,12]))
  fid[9,3] <- (1-Fid)*(move.primary[9,3]/(move.primary[9,1]+move.primary[9,2]+move.primary[9,3]+move.primary[9,4]+move.primary[9,5]+move.primary[9,6]+move.primary[9,7]+move.primary[9,8]+move.primary[9,10]+move.primary[9,11]+move.primary[9,12]))
  fid[9,4] <- (1-Fid)*(move.primary[9,4]/(move.primary[9,1]+move.primary[9,2]+move.primary[9,3]+move.primary[9,4]+move.primary[9,5]+move.primary[9,6]+move.primary[9,7]+move.primary[9,8]+move.primary[9,10]+move.primary[9,11]+move.primary[9,12]))
  fid[9,5] <- (1-Fid)*(move.primary[9,5]/(move.primary[9,1]+move.primary[9,2]+move.primary[9,3]+move.primary[9,4]+move.primary[9,5]+move.primary[9,6]+move.primary[9,7]+move.primary[9,8]+move.primary[9,10]+move.primary[9,11]+move.primary[9,12]))
  fid[9,6] <- (1-Fid)*(move.primary[9,6]/(move.primary[9,1]+move.primary[9,2]+move.primary[9,3]+move.primary[9,4]+move.primary[9,5]+move.primary[9,6]+move.primary[9,7]+move.primary[9,8]+move.primary[9,10]+move.primary[9,11]+move.primary[9,12]))
  fid[9,7] <- (1-Fid)*(move.primary[9,7]/(move.primary[9,1]+move.primary[9,2]+move.primary[9,3]+move.primary[9,4]+move.primary[9,5]+move.primary[9,6]+move.primary[9,7]+move.primary[9,8]+move.primary[9,10]+move.primary[9,11]+move.primary[9,12]))
  fid[9,8] <- (1-Fid)*(move.primary[9,8]/(move.primary[9,1]+move.primary[9,2]+move.primary[9,3]+move.primary[9,4]+move.primary[9,5]+move.primary[9,6]+move.primary[9,7]+move.primary[9,8]+move.primary[9,10]+move.primary[9,11]+move.primary[9,12]))
  fid[9,10] <- (1-Fid)*(move.primary[9,10]/(move.primary[9,1]+move.primary[9,2]+move.primary[9,3]+move.primary[9,4]+move.primary[9,5]+move.primary[9,6]+move.primary[9,7]+move.primary[9,8]+move.primary[9,10]+move.primary[9,11]+move.primary[9,12]))
  fid[9,11] <- (1-Fid)*(move.primary[9,11]/(move.primary[9,1]+move.primary[9,2]+move.primary[9,3]+move.primary[9,4]+move.primary[9,5]+move.primary[9,6]+move.primary[9,7]+move.primary[9,8]+move.primary[9,10]+move.primary[9,11]+move.primary[9,12]))
  fid[9,12] <- (1-Fid)*(move.primary[9,12]/(move.primary[9,1]+move.primary[9,2]+move.primary[9,3]+move.primary[9,4]+move.primary[9,5]+move.primary[9,6]+move.primary[9,7]+move.primary[9,8]+move.primary[9,10]+move.primary[9,11]+move.primary[9,12]))
  fid[9,9] <- Fid
  
  fid[10,1] <- (1-Fid)*(move.primary[10,1]/(move.primary[10,1]+move.primary[10,2]+move.primary[10,3]+move.primary[10,4]+move.primary[10,5]+move.primary[10,6]+move.primary[10,7]+move.primary[10,8]+move.primary[10,9]+move.primary[10,11]+move.primary[10,12]))
  fid[10,2] <- (1-Fid)*(move.primary[10,2]/(move.primary[10,1]+move.primary[10,2]+move.primary[10,3]+move.primary[10,4]+move.primary[10,5]+move.primary[10,6]+move.primary[10,7]+move.primary[10,8]+move.primary[10,9]+move.primary[10,11]+move.primary[10,12]))
  fid[10,3] <- (1-Fid)*(move.primary[10,3]/(move.primary[10,1]+move.primary[10,2]+move.primary[10,3]+move.primary[10,4]+move.primary[10,5]+move.primary[10,6]+move.primary[10,7]+move.primary[10,8]+move.primary[10,9]+move.primary[10,11]+move.primary[10,12]))
  fid[10,4] <- (1-Fid)*(move.primary[10,4]/(move.primary[10,1]+move.primary[10,2]+move.primary[10,3]+move.primary[10,4]+move.primary[10,5]+move.primary[10,6]+move.primary[10,7]+move.primary[10,8]+move.primary[10,9]+move.primary[10,11]+move.primary[10,12]))
  fid[10,5] <- (1-Fid)*(move.primary[10,5]/(move.primary[10,1]+move.primary[10,2]+move.primary[10,3]+move.primary[10,4]+move.primary[10,5]+move.primary[10,6]+move.primary[10,7]+move.primary[10,8]+move.primary[10,9]+move.primary[10,11]+move.primary[10,12]))
  fid[10,6] <- (1-Fid)*(move.primary[10,6]/(move.primary[10,1]+move.primary[10,2]+move.primary[10,3]+move.primary[10,4]+move.primary[10,5]+move.primary[10,6]+move.primary[10,7]+move.primary[10,8]+move.primary[10,9]+move.primary[10,11]+move.primary[10,12]))
  fid[10,7] <- (1-Fid)*(move.primary[10,7]/(move.primary[10,1]+move.primary[10,2]+move.primary[10,3]+move.primary[10,4]+move.primary[10,5]+move.primary[10,6]+move.primary[10,7]+move.primary[10,8]+move.primary[10,9]+move.primary[10,11]+move.primary[10,12]))
  fid[10,8] <- (1-Fid)*(move.primary[10,8]/(move.primary[10,1]+move.primary[10,2]+move.primary[10,3]+move.primary[10,4]+move.primary[10,5]+move.primary[10,6]+move.primary[10,7]+move.primary[10,8]+move.primary[10,9]+move.primary[10,11]+move.primary[10,12]))
  fid[10,9] <- (1-Fid)*(move.primary[10,9]/(move.primary[10,1]+move.primary[10,2]+move.primary[10,3]+move.primary[10,4]+move.primary[10,5]+move.primary[10,6]+move.primary[10,7]+move.primary[10,8]+move.primary[10,9]+move.primary[10,11]+move.primary[10,12]))
  fid[10,11] <- (1-Fid)*(move.primary[10,11]/(move.primary[10,1]+move.primary[10,2]+move.primary[10,3]+move.primary[10,4]+move.primary[10,5]+move.primary[10,6]+move.primary[10,7]+move.primary[10,8]+move.primary[10,9]+move.primary[10,11]+move.primary[10,12]))
  fid[10,12] <- (1-Fid)*(move.primary[10,12]/(move.primary[10,1]+move.primary[10,2]+move.primary[10,3]+move.primary[10,4]+move.primary[10,5]+move.primary[10,6]+move.primary[10,7]+move.primary[10,8]+move.primary[10,9]+move.primary[10,11]+move.primary[10,12]))
  fid[10,10] <- Fid
  
  fid[11,1] <- (1-Fid)*(move.primary[11,1]/(move.primary[11,1]+move.primary[11,2]+move.primary[11,3]+move.primary[11,4]+move.primary[11,5]+move.primary[11,6]+move.primary[11,7]+move.primary[11,8]+move.primary[11,9]+move.primary[11,10]+move.primary[11,12]))
  fid[11,2] <- (1-Fid)*(move.primary[11,2]/(move.primary[11,1]+move.primary[11,2]+move.primary[11,3]+move.primary[11,4]+move.primary[11,5]+move.primary[11,6]+move.primary[11,7]+move.primary[11,8]+move.primary[11,9]+move.primary[11,10]+move.primary[11,12]))
  fid[11,3] <- (1-Fid)*(move.primary[11,3]/(move.primary[11,1]+move.primary[11,2]+move.primary[11,3]+move.primary[11,4]+move.primary[11,5]+move.primary[11,6]+move.primary[11,7]+move.primary[11,8]+move.primary[11,9]+move.primary[11,10]+move.primary[11,12]))
  fid[11,4] <- (1-Fid)*(move.primary[11,4]/(move.primary[11,1]+move.primary[11,2]+move.primary[11,3]+move.primary[11,4]+move.primary[11,5]+move.primary[11,6]+move.primary[11,7]+move.primary[11,8]+move.primary[11,9]+move.primary[11,10]+move.primary[11,12]))
  fid[11,5] <- (1-Fid)*(move.primary[11,5]/(move.primary[11,1]+move.primary[11,2]+move.primary[11,3]+move.primary[11,4]+move.primary[11,5]+move.primary[11,6]+move.primary[11,7]+move.primary[11,8]+move.primary[11,9]+move.primary[11,10]+move.primary[11,12]))
  fid[11,6] <- (1-Fid)*(move.primary[11,6]/(move.primary[11,1]+move.primary[11,2]+move.primary[11,3]+move.primary[11,4]+move.primary[11,5]+move.primary[11,6]+move.primary[11,7]+move.primary[11,8]+move.primary[11,9]+move.primary[11,10]+move.primary[11,12]))
  fid[11,7] <- (1-Fid)*(move.primary[11,7]/(move.primary[11,1]+move.primary[11,2]+move.primary[11,3]+move.primary[11,4]+move.primary[11,5]+move.primary[11,6]+move.primary[11,7]+move.primary[11,8]+move.primary[11,9]+move.primary[11,10]+move.primary[11,12]))
  fid[11,8] <- (1-Fid)*(move.primary[11,8]/(move.primary[11,1]+move.primary[11,2]+move.primary[11,3]+move.primary[11,4]+move.primary[11,5]+move.primary[11,6]+move.primary[11,7]+move.primary[11,8]+move.primary[11,9]+move.primary[11,10]+move.primary[11,12]))
  fid[11,9] <- (1-Fid)*(move.primary[11,9]/(move.primary[11,1]+move.primary[11,2]+move.primary[11,3]+move.primary[11,4]+move.primary[11,5]+move.primary[11,6]+move.primary[11,7]+move.primary[11,8]+move.primary[11,9]+move.primary[11,10]+move.primary[11,12]))
  fid[11,10] <- (1-Fid)*(move.primary[11,10]/(move.primary[11,1]+move.primary[11,2]+move.primary[11,3]+move.primary[11,4]+move.primary[11,5]+move.primary[11,6]+move.primary[11,7]+move.primary[11,8]+move.primary[11,9]+move.primary[11,10]+move.primary[11,12]))
  fid[11,12] <- (1-Fid)*(move.primary[11,12]/(move.primary[11,1]+move.primary[11,2]+move.primary[11,3]+move.primary[11,4]+move.primary[11,5]+move.primary[11,6]+move.primary[11,7]+move.primary[11,8]+move.primary[11,9]+move.primary[11,10]+move.primary[11,12]))
  fid[11,11] <- Fid
  
  fid[12,1] <- (1-Fid)*(move.primary[12,1]/(move.primary[12,1]+move.primary[12,2]+move.primary[12,3]+move.primary[12,4]+move.primary[12,5]+move.primary[12,6]+move.primary[12,7]+move.primary[12,8]+move.primary[12,9]+move.primary[12,10]+move.primary[12,11]))
  fid[12,2] <- (1-Fid)*(move.primary[12,2]/(move.primary[12,1]+move.primary[12,2]+move.primary[12,3]+move.primary[12,4]+move.primary[12,5]+move.primary[12,6]+move.primary[12,7]+move.primary[12,8]+move.primary[12,9]+move.primary[12,10]+move.primary[12,11]))
  fid[12,3] <- (1-Fid)*(move.primary[12,3]/(move.primary[12,1]+move.primary[12,2]+move.primary[12,3]+move.primary[12,4]+move.primary[12,5]+move.primary[12,6]+move.primary[12,7]+move.primary[12,8]+move.primary[12,9]+move.primary[12,10]+move.primary[12,11]))
  fid[12,4] <- (1-Fid)*(move.primary[12,4]/(move.primary[12,1]+move.primary[12,2]+move.primary[12,3]+move.primary[12,4]+move.primary[12,5]+move.primary[12,6]+move.primary[12,7]+move.primary[12,8]+move.primary[12,9]+move.primary[12,10]+move.primary[12,11]))
  fid[12,5] <- (1-Fid)*(move.primary[12,5]/(move.primary[12,1]+move.primary[12,2]+move.primary[12,3]+move.primary[12,4]+move.primary[12,5]+move.primary[12,6]+move.primary[12,7]+move.primary[12,8]+move.primary[12,9]+move.primary[12,10]+move.primary[12,11]))
  fid[12,6] <- (1-Fid)*(move.primary[12,6]/(move.primary[12,1]+move.primary[12,2]+move.primary[12,3]+move.primary[12,4]+move.primary[12,5]+move.primary[12,6]+move.primary[12,7]+move.primary[12,8]+move.primary[12,9]+move.primary[12,10]+move.primary[12,11]))
  fid[12,7] <- (1-Fid)*(move.primary[12,7]/(move.primary[12,1]+move.primary[12,2]+move.primary[12,3]+move.primary[12,4]+move.primary[12,5]+move.primary[12,6]+move.primary[12,7]+move.primary[12,8]+move.primary[12,9]+move.primary[12,10]+move.primary[12,11]))
  fid[12,8] <- (1-Fid)*(move.primary[12,8]/(move.primary[12,1]+move.primary[12,2]+move.primary[12,3]+move.primary[12,4]+move.primary[12,5]+move.primary[12,6]+move.primary[12,7]+move.primary[12,8]+move.primary[12,9]+move.primary[12,10]+move.primary[12,11]))
  fid[12,9] <- (1-Fid)*(move.primary[12,9]/(move.primary[12,1]+move.primary[12,2]+move.primary[12,3]+move.primary[12,4]+move.primary[12,5]+move.primary[12,6]+move.primary[12,7]+move.primary[12,8]+move.primary[12,9]+move.primary[12,10]+move.primary[12,11]))
  fid[12,10] <- (1-Fid)*(move.primary[12,10]/(move.primary[12,1]+move.primary[12,2]+move.primary[12,3]+move.primary[12,4]+move.primary[12,5]+move.primary[12,6]+move.primary[12,7]+move.primary[12,8]+move.primary[12,9]+move.primary[12,10]+move.primary[12,11]))
  fid[12,11] <- (1-Fid)*(move.primary[12,11]/(move.primary[12,1]+move.primary[12,2]+move.primary[12,3]+move.primary[12,4]+move.primary[12,5]+move.primary[12,6]+move.primary[12,7]+move.primary[12,8]+move.primary[12,9]+move.primary[12,10]+move.primary[12,11]))
  fid[12,12] <- Fid
  
  # transition probability matrices between periods
  omega1 <- 1/(1+exp(-param[41]))
  omega2 <- 1/(1+exp(-param[42]))
  
  gamma <- array(0, dim=c(26, 26, T-1))
  gamma[1,1,] <- 1-rstar[2:T] #rstar is recruitment into the breeding population - all ponds are the same
  gamma[1,2,] <- rstar[2:T]*alpha[1]
  gamma[1,3,] <- rstar[2:T]*alpha[2]
  gamma[1,4,] <- rstar[2:T]*alpha[3]
  gamma[1,5,] <- rstar[2:T]*alpha[4]
  gamma[1,6,] <- rstar[2:T]*alpha[5]
  gamma[1,7,] <- rstar[2:T]*alpha[6]
  gamma[1,8,] <- rstar[2:T]*alpha[7]
  gamma[1,9,] <- rstar[2:T]*alpha[8]
  gamma[1,10,] <- rstar[2:T]*alpha[9]
  gamma[1,11,] <- rstar[2:T]*alpha[10]
  gamma[1,12,] <- rstar[2:T]*alpha[11]
  gamma[1,13,] <- rstar[2:T]*alpha[12]
  
  
  gamma[2,2,] <- s * fid[1,1] * (1-omega2)
  gamma[2,3,] <- s * fid[1,2] * (1-omega2)
  gamma[2,4,] <- s * fid[1,3] * (1-omega2)
  gamma[2,5,] <- s * fid[1,4] * (1-omega2)
  gamma[2,6,] <- s * fid[1,5] * (1-omega2)
  gamma[2,7,] <- s * fid[1,6] * (1-omega2)
  gamma[2,8,] <- s * fid[1,7] * (1-omega2)
  gamma[2,9,] <- s * fid[1,8] * (1-omega2)
  gamma[2,10,] <- s * fid[1,9] * (1-omega2)
  gamma[2,11,] <- s * fid[1,10] * (1-omega2)
  gamma[2,12,] <- s * fid[1,11] * (1-omega2)
  gamma[2,13,] <- s * fid[1,12] * (1-omega2)
  gamma[2,14,] <- s * fid[1,1] * omega2
  gamma[2,15,] <- s * fid[1,2] * omega2
  gamma[2,16,] <- s * fid[1,3] * omega2
  gamma[2,17,] <- s * fid[1,4] * omega2
  gamma[2,18,] <- s * fid[1,5] * omega2
  gamma[2,19,] <- s * fid[1,6] * omega2
  gamma[2,20,] <- s * fid[1,7] * omega2
  gamma[2,21,] <- s * fid[1,8] * omega2
  gamma[2,22,] <- s * fid[1,9] * omega2
  gamma[2,23,] <- s * fid[1,10] * omega2
  gamma[2,24,] <- s * fid[1,11] * omega2
  gamma[2,25,] <- s * fid[1,12] * omega2
  gamma[2,26,] <- 1-s
  
  gamma[3,2,] <- s * fid[2,1] * (1-omega2)
  gamma[3,3,] <- s * fid[2,2] * (1-omega2)
  gamma[3,4,] <- s * fid[2,3] * (1-omega2)
  gamma[3,5,] <- s * fid[2,4] * (1-omega2)
  gamma[3,6,] <- s * fid[2,5] * (1-omega2)
  gamma[3,7,] <- s * fid[2,6] * (1-omega2)
  gamma[3,8,] <- s * fid[2,7] * (1-omega2)
  gamma[3,9,] <- s * fid[2,8] * (1-omega2)
  gamma[3,10,] <- s * fid[2,9] * (1-omega2)
  gamma[3,11,] <- s * fid[2,10] * (1-omega2)
  gamma[3,12,] <- s * fid[2,11] * (1-omega2)
  gamma[3,13,] <- s * fid[2,12] * (1-omega2)
  gamma[3,14,] <- s * fid[2,1] * omega2
  gamma[3,15,] <- s * fid[2,2] * omega2
  gamma[3,16,] <- s * fid[2,3] * omega2
  gamma[3,17,] <- s * fid[2,4] * omega2
  gamma[3,18,] <- s * fid[2,5] * omega2
  gamma[3,19,] <- s * fid[2,6] * omega2
  gamma[3,20,] <- s * fid[2,7] * omega2
  gamma[3,21,] <- s * fid[2,8] * omega2
  gamma[3,22,] <- s * fid[2,9] * omega2
  gamma[3,23,] <- s * fid[2,10] * omega2
  gamma[3,24,] <- s * fid[2,11] * omega2
  gamma[3,25,] <- s * fid[2,12] * omega2
  gamma[3,26,] <- 1-s
  
  gamma[4,2,] <- s * fid[3,1] * (1-omega2)
  gamma[4,3,] <- s * fid[3,2] * (1-omega2)
  gamma[4,4,] <- s * fid[3,3] * (1-omega2)
  gamma[4,5,] <- s * fid[3,4] * (1-omega2)
  gamma[4,6,] <- s * fid[3,5] * (1-omega2)
  gamma[4,7,] <- s * fid[3,6] * (1-omega2)
  gamma[4,8,] <- s * fid[3,7] * (1-omega2)
  gamma[4,9,] <- s * fid[3,8] * (1-omega2)
  gamma[4,10,] <- s * fid[3,9] * (1-omega2)
  gamma[4,11,] <- s * fid[3,10] * (1-omega2)
  gamma[4,12,] <- s * fid[3,11] * (1-omega2)
  gamma[4,13,] <- s * fid[3,12] * (1-omega2)
  gamma[4,14,] <- s * fid[3,1] * omega2
  gamma[4,15,] <- s * fid[3,2] * omega2
  gamma[4,16,] <- s * fid[3,3] * omega2
  gamma[4,17,] <- s * fid[3,4] * omega2
  gamma[4,18,] <- s * fid[3,5] * omega2
  gamma[4,19,] <- s * fid[3,6] * omega2
  gamma[4,20,] <- s * fid[3,7] * omega2
  gamma[4,21,] <- s * fid[3,8] * omega2
  gamma[4,22,] <- s * fid[3,9] * omega2
  gamma[4,23,] <- s * fid[3,10] * omega2
  gamma[4,24,] <- s * fid[3,11] * omega2
  gamma[4,25,] <- s * fid[3,12] * omega2
  gamma[4,26,] <- 1-s
  
  gamma[5,2,] <- s * fid[4,1] * (1-omega2)
  gamma[5,3,] <- s * fid[4,2] * (1-omega2)
  gamma[5,4,] <- s * fid[4,3] * (1-omega2)
  gamma[5,5,] <- s * fid[4,4] * (1-omega2)
  gamma[5,6,] <- s * fid[4,5] * (1-omega2)
  gamma[5,7,] <- s * fid[4,6] * (1-omega2)
  gamma[5,8,] <- s * fid[4,7] * (1-omega2)
  gamma[5,9,] <- s * fid[4,8] * (1-omega2)
  gamma[5,10,] <- s * fid[4,9] * (1-omega2)
  gamma[5,11,] <- s * fid[4,10] * (1-omega2)
  gamma[5,12,] <- s * fid[4,11] * (1-omega2)
  gamma[5,13,] <- s * fid[4,12] * (1-omega2)
  gamma[5,14,] <- s * fid[4,1] * omega2
  gamma[5,15,] <- s * fid[4,2] * omega2
  gamma[5,16,] <- s * fid[4,3] * omega2
  gamma[5,17,] <- s * fid[4,4] * omega2
  gamma[5,18,] <- s * fid[4,5] * omega2
  gamma[5,19,] <- s * fid[4,6] * omega2
  gamma[5,20,] <- s * fid[4,7] * omega2
  gamma[5,21,] <- s * fid[4,8] * omega2
  gamma[5,22,] <- s * fid[4,9] * omega2
  gamma[5,23,] <- s * fid[4,10] * omega2
  gamma[5,24,] <- s * fid[4,11] * omega2
  gamma[5,25,] <- s * fid[4,12] * omega2
  gamma[5,26,] <- 1-s
  
  gamma[6,2,] <- s * fid[5,1] * (1-omega2)
  gamma[6,3,] <- s * fid[5,2] * (1-omega2)
  gamma[6,4,] <- s * fid[5,3] * (1-omega2)
  gamma[6,5,] <- s * fid[5,4] * (1-omega2)
  gamma[6,6,] <- s * fid[5,5] * (1-omega2)
  gamma[6,7,] <- s * fid[5,6] * (1-omega2)
  gamma[6,8,] <- s * fid[5,7] * (1-omega2)
  gamma[6,9,] <- s * fid[5,8] * (1-omega2)
  gamma[6,10,] <- s * fid[5,9] * (1-omega2)
  gamma[6,11,] <- s * fid[5,10] * (1-omega2)
  gamma[6,12,] <- s * fid[5,11] * (1-omega2)
  gamma[6,13,] <- s * fid[5,12] * (1-omega2)
  gamma[6,14,] <- s * fid[5,1] * omega2
  gamma[6,15,] <- s * fid[5,2] * omega2
  gamma[6,16,] <- s * fid[5,3] * omega2
  gamma[6,17,] <- s * fid[5,4] * omega2
  gamma[6,18,] <- s * fid[5,5] * omega2
  gamma[6,19,] <- s * fid[5,6] * omega2
  gamma[6,20,] <- s * fid[5,7] * omega2
  gamma[6,21,] <- s * fid[5,8] * omega2
  gamma[6,22,] <- s * fid[5,9] * omega2
  gamma[6,23,] <- s * fid[5,10] * omega2
  gamma[6,24,] <- s * fid[5,11] * omega2
  gamma[6,25,] <- s * fid[5,12] * omega2
  gamma[6,26,] <- 1-s
  
  gamma[7,2,] <- s * fid[6,1] * (1-omega2)
  gamma[7,3,] <- s * fid[6,2] * (1-omega2)
  gamma[7,4,] <- s * fid[6,3] * (1-omega2)
  gamma[7,5,] <- s * fid[6,4] * (1-omega2)
  gamma[7,6,] <- s * fid[6,5] * (1-omega2)
  gamma[7,7,] <- s * fid[6,6] * (1-omega2)
  gamma[7,8,] <- s * fid[6,7] * (1-omega2)
  gamma[7,9,] <- s * fid[6,8] * (1-omega2)
  gamma[7,10,] <- s * fid[6,9] * (1-omega2)
  gamma[7,11,] <- s * fid[6,10] * (1-omega2)
  gamma[7,12,] <- s * fid[6,11] * (1-omega2)
  gamma[7,13,] <- s * fid[6,12] * (1-omega2)
  gamma[7,14,] <- s * fid[6,1] * omega2
  gamma[7,15,] <- s * fid[6,2] * omega2
  gamma[7,16,] <- s * fid[6,3] * omega2
  gamma[7,17,] <- s * fid[6,4] * omega2
  gamma[7,18,] <- s * fid[6,5] * omega2
  gamma[7,19,] <- s * fid[6,6] * omega2
  gamma[7,20,] <- s * fid[6,7] * omega2
  gamma[7,21,] <- s * fid[6,8] * omega2
  gamma[7,22,] <- s * fid[6,9] * omega2
  gamma[7,23,] <- s * fid[6,10] * omega2
  gamma[7,24,] <- s * fid[6,11] * omega2
  gamma[7,25,] <- s * fid[6,12] * omega2
  gamma[7,26,] <- 1-s
  
  gamma[8,2,] <- s * fid[7,1] * (1-omega2)
  gamma[8,3,] <- s * fid[7,2] * (1-omega2)
  gamma[8,4,] <- s * fid[7,3] * (1-omega2)
  gamma[8,5,] <- s * fid[7,4] * (1-omega2)
  gamma[8,6,] <- s * fid[7,5] * (1-omega2)
  gamma[8,7,] <- s * fid[7,6] * (1-omega2)
  gamma[8,8,] <- s * fid[7,7] * (1-omega2)
  gamma[8,9,] <- s * fid[7,8] * (1-omega2)
  gamma[8,10,] <- s * fid[7,9] * (1-omega2)
  gamma[8,11,] <- s * fid[7,10] * (1-omega2)
  gamma[8,12,] <- s * fid[7,11] * (1-omega2)
  gamma[8,13,] <- s * fid[7,12] * (1-omega2)
  gamma[8,14,] <- s * fid[7,1] * omega2
  gamma[8,15,] <- s * fid[7,2] * omega2
  gamma[8,16,] <- s * fid[7,3] * omega2
  gamma[8,17,] <- s * fid[7,4] * omega2
  gamma[8,18,] <- s * fid[7,5] * omega2
  gamma[8,19,] <- s * fid[7,6] * omega2
  gamma[8,20,] <- s * fid[7,7] * omega2
  gamma[8,21,] <- s * fid[7,8] * omega2
  gamma[8,22,] <- s * fid[7,9] * omega2
  gamma[8,23,] <- s * fid[7,10] * omega2
  gamma[8,24,] <- s * fid[7,11] * omega2
  gamma[8,25,] <- s * fid[7,12] * omega2
  gamma[8,26,] <- 1-s
  
  gamma[9,2,] <- s * fid[8,1] * (1-omega2)
  gamma[9,3,] <- s * fid[8,2] * (1-omega2)
  gamma[9,4,] <- s * fid[8,3] * (1-omega2)
  gamma[9,5,] <- s * fid[8,4] * (1-omega2)
  gamma[9,6,] <- s * fid[8,5] * (1-omega2)
  gamma[9,7,] <- s * fid[8,6] * (1-omega2)
  gamma[9,8,] <- s * fid[8,7] * (1-omega2)
  gamma[9,9,] <- s * fid[8,8] * (1-omega2)
  gamma[9,10,] <- s * fid[8,9] * (1-omega2)
  gamma[9,11,] <- s * fid[8,10] * (1-omega2)
  gamma[9,12,] <- s * fid[8,11] * (1-omega2)
  gamma[9,13,] <- s * fid[8,12] * (1-omega2)
  gamma[9,14,] <- s * fid[8,1] * omega2
  gamma[9,15,] <- s * fid[8,2] * omega2
  gamma[9,16,] <- s * fid[8,3] * omega2
  gamma[9,17,] <- s * fid[8,4] * omega2
  gamma[9,18,] <- s * fid[8,5] * omega2
  gamma[9,19,] <- s * fid[8,6] * omega2
  gamma[9,20,] <- s * fid[8,7] * omega2
  gamma[9,21,] <- s * fid[8,8] * omega2
  gamma[9,22,] <- s * fid[8,9] * omega2
  gamma[9,23,] <- s * fid[8,10] * omega2
  gamma[9,24,] <- s * fid[8,11] * omega2
  gamma[9,25,] <- s * fid[8,12] * omega2
  gamma[9,26,] <- 1-s
  
  gamma[10,2,] <- s * fid[9,1] * (1-omega2)
  gamma[10,3,] <- s * fid[9,2] * (1-omega2)
  gamma[10,4,] <- s * fid[9,3] * (1-omega2)
  gamma[10,5,] <- s * fid[9,4] * (1-omega2)
  gamma[10,6,] <- s * fid[9,5] * (1-omega2)
  gamma[10,7,] <- s * fid[9,6] * (1-omega2)
  gamma[10,8,] <- s * fid[9,7] * (1-omega2)
  gamma[10,9,] <- s * fid[9,8] * (1-omega2)
  gamma[10,10,] <- s * fid[9,9] * (1-omega2)
  gamma[10,11,] <- s * fid[9,10] * (1-omega2)
  gamma[10,12,] <- s * fid[9,11] * (1-omega2)
  gamma[10,13,] <- s * fid[9,12] * (1-omega2)
  gamma[10,14,] <- s * fid[9,1] * omega2
  gamma[10,15,] <- s * fid[9,2] * omega2
  gamma[10,16,] <- s * fid[9,3] * omega2
  gamma[10,17,] <- s * fid[9,4] * omega2
  gamma[10,18,] <- s * fid[9,5] * omega2
  gamma[10,19,] <- s * fid[9,6] * omega2
  gamma[10,20,] <- s * fid[9,7] * omega2
  gamma[10,21,] <- s * fid[9,8] * omega2
  gamma[10,22,] <- s * fid[9,9] * omega2
  gamma[10,23,] <- s * fid[9,10] * omega2
  gamma[10,24,] <- s * fid[9,11] * omega2
  gamma[10,25,] <- s * fid[9,12] * omega2
  gamma[10,26,] <- 1-s
  
  gamma[11,2,] <- s * fid[10,1] * (1-omega2)
  gamma[11,3,] <- s * fid[10,2] * (1-omega2)
  gamma[11,4,] <- s * fid[10,3] * (1-omega2)
  gamma[11,5,] <- s * fid[10,4] * (1-omega2)
  gamma[11,6,] <- s * fid[10,5] * (1-omega2)
  gamma[11,7,] <- s * fid[10,6] * (1-omega2)
  gamma[11,8,] <- s * fid[10,7] * (1-omega2)
  gamma[11,9,] <- s * fid[10,8] * (1-omega2)
  gamma[11,10,] <- s * fid[10,9] * (1-omega2)
  gamma[11,11,] <- s * fid[10,10] * (1-omega2)
  gamma[11,12,] <- s * fid[10,11] * (1-omega2)
  gamma[11,13,] <- s * fid[10,12] * (1-omega2)
  gamma[11,14,] <- s * fid[10,1] * omega2
  gamma[11,15,] <- s * fid[10,2] * omega2
  gamma[11,16,] <- s * fid[10,3] * omega2
  gamma[11,17,] <- s * fid[10,4] * omega2
  gamma[11,18,] <- s * fid[10,5] * omega2
  gamma[11,19,] <- s * fid[10,6] * omega2
  gamma[11,20,] <- s * fid[10,7] * omega2
  gamma[11,21,] <- s * fid[10,8] * omega2
  gamma[11,22,] <- s * fid[10,9] * omega2
  gamma[11,23,] <- s * fid[10,10] * omega2
  gamma[11,24,] <- s * fid[10,11] * omega2
  gamma[11,25,] <- s * fid[10,12] * omega2
  gamma[11,26,] <- 1-s
  
  gamma[12,2,] <- s * fid[11,1] * (1-omega2)
  gamma[12,3,] <- s * fid[11,2] * (1-omega2)
  gamma[12,4,] <- s * fid[11,3] * (1-omega2)
  gamma[12,5,] <- s * fid[11,4] * (1-omega2)
  gamma[12,6,] <- s * fid[11,5] * (1-omega2)
  gamma[12,7,] <- s * fid[11,6] * (1-omega2)
  gamma[12,8,] <- s * fid[11,7] * (1-omega2)
  gamma[12,9,] <- s * fid[11,8] * (1-omega2)
  gamma[12,10,] <- s * fid[11,9] * (1-omega2)
  gamma[12,11,] <- s * fid[11,10] * (1-omega2)
  gamma[12,12,] <- s * fid[11,11] * (1-omega2)
  gamma[12,13,] <- s * fid[11,12] * (1-omega2)
  gamma[12,14,] <- s * fid[11,1] * omega2
  gamma[12,15,] <- s * fid[11,2] * omega2
  gamma[12,16,] <- s * fid[11,3] * omega2
  gamma[12,17,] <- s * fid[11,4] * omega2
  gamma[12,18,] <- s * fid[11,5] * omega2
  gamma[12,19,] <- s * fid[11,6] * omega2
  gamma[12,20,] <- s * fid[11,7] * omega2
  gamma[12,21,] <- s * fid[11,8] * omega2
  gamma[12,22,] <- s * fid[11,9] * omega2
  gamma[12,23,] <- s * fid[11,10] * omega2
  gamma[12,24,] <- s * fid[11,11] * omega2
  gamma[12,25,] <- s * fid[11,12] * omega2
  gamma[12,26,] <- 1-s
  
  gamma[13,2,] <- s * fid[12,1] * (1-omega2)
  gamma[13,3,] <- s * fid[12,2] * (1-omega2)
  gamma[13,4,] <- s * fid[12,3] * (1-omega2)
  gamma[13,5,] <- s * fid[12,4] * (1-omega2)
  gamma[13,6,] <- s * fid[12,5] * (1-omega2)
  gamma[13,7,] <- s * fid[12,6] * (1-omega2)
  gamma[13,8,] <- s * fid[12,7] * (1-omega2)
  gamma[13,9,] <- s * fid[12,8] * (1-omega2)
  gamma[13,10,] <- s * fid[12,9] * (1-omega2)
  gamma[13,11,] <- s * fid[12,10] * (1-omega2)
  gamma[13,12,] <- s * fid[12,11] * (1-omega2)
  gamma[13,13,] <- s * fid[12,12] * (1-omega2)
  gamma[13,14,] <- s * fid[12,1] * omega2
  gamma[13,15,] <- s * fid[12,2] * omega2
  gamma[13,16,] <- s * fid[12,3] * omega2
  gamma[13,17,] <- s * fid[12,4] * omega2
  gamma[13,18,] <- s * fid[12,5] * omega2
  gamma[13,19,] <- s * fid[12,6] * omega2
  gamma[13,20,] <- s * fid[12,7] * omega2
  gamma[13,21,] <- s * fid[12,8] * omega2
  gamma[13,22,] <- s * fid[12,9] * omega2
  gamma[13,23,] <- s * fid[12,10] * omega2
  gamma[13,24,] <- s * fid[12,11] * omega2
  gamma[13,25,] <- s * fid[12,12] * omega2
  gamma[13,26,] <- 1-s
  
  gamma[14,2,] <- s * fid[1,1] * (1-omega1)
  gamma[14,3,] <- s * fid[1,2] * (1-omega1)
  gamma[14,4,] <- s * fid[1,3] * (1-omega1)
  gamma[14,5,] <- s * fid[1,4] * (1-omega1)
  gamma[14,6,] <- s * fid[1,5] * (1-omega1)
  gamma[14,7,] <- s * fid[1,6] * (1-omega1)
  gamma[14,8,] <- s * fid[1,7] * (1-omega1)
  gamma[14,9,] <- s * fid[1,8] * (1-omega1)
  gamma[14,10,] <- s * fid[1,9] * (1-omega1)
  gamma[14,11,] <- s * fid[1,10] * (1-omega1)
  gamma[14,12,] <- s * fid[1,11] * (1-omega1)
  gamma[14,13,] <- s * fid[1,12] * (1-omega1)
  gamma[14,14,] <- s * fid[1,1] * omega1
  gamma[14,15,] <- s * fid[1,2] * omega1
  gamma[14,16,] <- s * fid[1,3] * omega1
  gamma[14,17,] <- s * fid[1,4] * omega1
  gamma[14,18,] <- s * fid[1,5] * omega1
  gamma[14,19,] <- s * fid[1,6] * omega1
  gamma[14,20,] <- s * fid[1,7] * omega1
  gamma[14,21,] <- s * fid[1,8] * omega1
  gamma[14,22,] <- s * fid[1,9] * omega1
  gamma[14,23,] <- s * fid[1,10] * omega1
  gamma[14,24,] <- s * fid[1,11] * omega1
  gamma[14,25,] <- s * fid[1,12] * omega1
  gamma[14,26,] <- 1-s
  
  gamma[15,2,] <- s * fid[2,1] * (1-omega1)
  gamma[15,3,] <- s * fid[2,2] * (1-omega1)
  gamma[15,4,] <- s * fid[2,3] * (1-omega1)
  gamma[15,5,] <- s * fid[2,4] * (1-omega1)
  gamma[15,6,] <- s * fid[2,5] * (1-omega1)
  gamma[15,7,] <- s * fid[2,6] * (1-omega1)
  gamma[15,8,] <- s * fid[2,7] * (1-omega1)
  gamma[15,9,] <- s * fid[2,8] * (1-omega1)
  gamma[15,10,] <- s * fid[2,9] * (1-omega1)
  gamma[15,11,] <- s * fid[2,10] * (1-omega1)
  gamma[15,12,] <- s * fid[2,11] * (1-omega1)
  gamma[15,13,] <- s * fid[2,12] * (1-omega1)
  gamma[15,14,] <- s * fid[2,1] * omega1
  gamma[15,15,] <- s * fid[2,2] * omega1
  gamma[15,16,] <- s * fid[2,3] * omega1
  gamma[15,17,] <- s * fid[2,4] * omega1
  gamma[15,18,] <- s * fid[2,5] * omega1
  gamma[15,19,] <- s * fid[2,6] * omega1
  gamma[15,20,] <- s * fid[2,7] * omega1
  gamma[15,21,] <- s * fid[2,8] * omega1
  gamma[15,22,] <- s * fid[2,9] * omega1
  gamma[15,23,] <- s * fid[2,10] * omega1
  gamma[15,24,] <- s * fid[2,11] * omega1
  gamma[15,25,] <- s * fid[2,12] * omega1
  gamma[15,26,] <- 1-s
  
  gamma[16,2,] <- s * fid[3,1] * (1-omega1)
  gamma[16,3,] <- s * fid[3,2] * (1-omega1)
  gamma[16,4,] <- s * fid[3,3] * (1-omega1)
  gamma[16,5,] <- s * fid[3,4] * (1-omega1)
  gamma[16,6,] <- s * fid[3,5] * (1-omega1)
  gamma[16,7,] <- s * fid[3,6] * (1-omega1)
  gamma[16,8,] <- s * fid[3,7] * (1-omega1)
  gamma[16,9,] <- s * fid[3,8] * (1-omega1)
  gamma[16,10,] <- s * fid[3,9] * (1-omega1)
  gamma[16,11,] <- s * fid[3,10] * (1-omega1)
  gamma[16,12,] <- s * fid[3,11] * (1-omega1)
  gamma[16,13,] <- s * fid[3,12] * (1-omega1)
  gamma[16,14,] <- s * fid[3,1] * omega1
  gamma[16,15,] <- s * fid[3,2] * omega1
  gamma[16,16,] <- s * fid[3,3] * omega1
  gamma[16,17,] <- s * fid[3,4] * omega1
  gamma[16,18,] <- s * fid[3,5] * omega1
  gamma[16,19,] <- s * fid[3,6] * omega1
  gamma[16,20,] <- s * fid[3,7] * omega1
  gamma[16,21,] <- s * fid[3,8] * omega1
  gamma[16,22,] <- s * fid[3,9] * omega1
  gamma[16,23,] <- s * fid[3,10] * omega1
  gamma[16,24,] <- s * fid[3,11] * omega1
  gamma[16,25,] <- s * fid[3,12] * omega1
  gamma[16,26,] <- 1-s
  
  gamma[17,2,] <- s * fid[4,1] * (1-omega1)
  gamma[17,3,] <- s * fid[4,2] * (1-omega1)
  gamma[17,4,] <- s * fid[4,3] * (1-omega1)
  gamma[17,5,] <- s * fid[4,4] * (1-omega1)
  gamma[17,6,] <- s * fid[4,5] * (1-omega1)
  gamma[17,7,] <- s * fid[4,6] * (1-omega1)
  gamma[17,8,] <- s * fid[4,7] * (1-omega1)
  gamma[17,9,] <- s * fid[4,8] * (1-omega1)
  gamma[17,10,] <- s * fid[4,9] * (1-omega1)
  gamma[17,11,] <- s * fid[4,10] * (1-omega1)
  gamma[17,12,] <- s * fid[4,11] * (1-omega1)
  gamma[17,13,] <- s * fid[4,12] * (1-omega1)
  gamma[17,14,] <- s * fid[4,1] * omega1
  gamma[17,15,] <- s * fid[4,2] * omega1
  gamma[17,16,] <- s * fid[4,3] * omega1
  gamma[17,17,] <- s * fid[4,4] * omega1
  gamma[17,18,] <- s * fid[4,5] * omega1
  gamma[17,19,] <- s * fid[4,6] * omega1
  gamma[17,20,] <- s * fid[4,7] * omega1
  gamma[17,21,] <- s * fid[4,8] * omega1
  gamma[17,22,] <- s * fid[4,9] * omega1
  gamma[17,23,] <- s * fid[4,10] * omega1
  gamma[17,24,] <- s * fid[4,11] * omega1
  gamma[17,25,] <- s * fid[4,12] * omega1
  gamma[17,26,] <- 1-s
  
  gamma[18,2,] <- s * fid[5,1] * (1-omega1)
  gamma[18,3,] <- s * fid[5,2] * (1-omega1)
  gamma[18,4,] <- s * fid[5,3] * (1-omega1)
  gamma[18,5,] <- s * fid[5,4] * (1-omega1)
  gamma[18,6,] <- s * fid[5,5] * (1-omega1)
  gamma[18,7,] <- s * fid[5,6] * (1-omega1)
  gamma[18,8,] <- s * fid[5,7] * (1-omega1)
  gamma[18,9,] <- s * fid[5,8] * (1-omega1)
  gamma[18,10,] <- s * fid[5,9] * (1-omega1)
  gamma[18,11,] <- s * fid[5,10] * (1-omega1)
  gamma[18,12,] <- s * fid[5,11] * (1-omega1)
  gamma[18,13,] <- s * fid[5,12] * (1-omega1)
  gamma[18,14,] <- s * fid[5,1] * omega1
  gamma[18,15,] <- s * fid[5,2] * omega1
  gamma[18,16,] <- s * fid[5,3] * omega1
  gamma[18,17,] <- s * fid[5,4] * omega1
  gamma[18,18,] <- s * fid[5,5] * omega1
  gamma[18,19,] <- s * fid[5,6] * omega1
  gamma[18,20,] <- s * fid[5,7] * omega1
  gamma[18,21,] <- s * fid[5,8] * omega1
  gamma[18,22,] <- s * fid[5,9] * omega1
  gamma[18,23,] <- s * fid[5,10] * omega1
  gamma[18,24,] <- s * fid[5,11] * omega1
  gamma[18,25,] <- s * fid[5,12] * omega1
  gamma[18,26,] <- 1-s
  
  gamma[19,2,] <- s * fid[6,1] * (1-omega1)
  gamma[19,3,] <- s * fid[6,2] * (1-omega1)
  gamma[19,4,] <- s * fid[6,3] * (1-omega1)
  gamma[19,5,] <- s * fid[6,4] * (1-omega1)
  gamma[19,6,] <- s * fid[6,5] * (1-omega1)
  gamma[19,7,] <- s * fid[6,6] * (1-omega1)
  gamma[19,8,] <- s * fid[6,7] * (1-omega1)
  gamma[19,9,] <- s * fid[6,8] * (1-omega1)
  gamma[19,10,] <- s * fid[6,9] * (1-omega1)
  gamma[19,11,] <- s * fid[6,10] * (1-omega1)
  gamma[19,12,] <- s * fid[6,11] * (1-omega1)
  gamma[19,13,] <- s * fid[6,12] * (1-omega1)
  gamma[19,14,] <- s * fid[6,1] * omega1
  gamma[19,15,] <- s * fid[6,2] * omega1
  gamma[19,16,] <- s * fid[6,3] * omega1
  gamma[19,17,] <- s * fid[6,4] * omega1
  gamma[19,18,] <- s * fid[6,5] * omega1
  gamma[19,19,] <- s * fid[6,6] * omega1
  gamma[19,20,] <- s * fid[6,7] * omega1
  gamma[19,21,] <- s * fid[6,8] * omega1
  gamma[19,22,] <- s * fid[6,9] * omega1
  gamma[19,23,] <- s * fid[6,10] * omega1
  gamma[19,24,] <- s * fid[6,11] * omega1
  gamma[19,25,] <- s * fid[6,12] * omega1
  gamma[19,26,] <- 1-s
  
  gamma[20,2,] <- s * fid[7,1] * (1-omega1)
  gamma[20,3,] <- s * fid[7,2] * (1-omega1)
  gamma[20,4,] <- s * fid[7,3] * (1-omega1)
  gamma[20,5,] <- s * fid[7,4] * (1-omega1)
  gamma[20,6,] <- s * fid[7,5] * (1-omega1)
  gamma[20,7,] <- s * fid[7,6] * (1-omega1)
  gamma[20,8,] <- s * fid[7,7] * (1-omega1)
  gamma[20,9,] <- s * fid[7,8] * (1-omega1)
  gamma[20,10,] <- s * fid[7,9] * (1-omega1)
  gamma[20,11,] <- s * fid[7,10] * (1-omega1)
  gamma[20,12,] <- s * fid[7,11] * (1-omega1)
  gamma[20,13,] <- s * fid[7,12] * (1-omega1)
  gamma[20,14,] <- s * fid[7,1] * omega1
  gamma[20,15,] <- s * fid[7,2] * omega1
  gamma[20,16,] <- s * fid[7,3] * omega1
  gamma[20,17,] <- s * fid[7,4] * omega1
  gamma[20,18,] <- s * fid[7,5] * omega1
  gamma[20,19,] <- s * fid[7,6] * omega1
  gamma[20,20,] <- s * fid[7,7] * omega1
  gamma[20,21,] <- s * fid[7,8] * omega1
  gamma[20,22,] <- s * fid[7,9] * omega1
  gamma[20,23,] <- s * fid[7,10] * omega1
  gamma[20,24,] <- s * fid[7,11] * omega1
  gamma[20,25,] <- s * fid[7,12] * omega1
  gamma[20,26,] <- 1-s
  
  gamma[21,2,] <- s * fid[8,1] * (1-omega1)
  gamma[21,3,] <- s * fid[8,2] * (1-omega1)
  gamma[21,4,] <- s * fid[8,3] * (1-omega1)
  gamma[21,5,] <- s * fid[8,4] * (1-omega1)
  gamma[21,6,] <- s * fid[8,5] * (1-omega1)
  gamma[21,7,] <- s * fid[8,6] * (1-omega1)
  gamma[21,8,] <- s * fid[8,7] * (1-omega1)
  gamma[21,9,] <- s * fid[8,8] * (1-omega1)
  gamma[21,10,] <- s * fid[8,9] * (1-omega1)
  gamma[21,11,] <- s * fid[8,10] * (1-omega1)
  gamma[21,12,] <- s * fid[8,11] * (1-omega1)
  gamma[21,13,] <- s * fid[8,12] * (1-omega1)
  gamma[21,14,] <- s * fid[8,1] * omega1
  gamma[21,15,] <- s * fid[8,2] * omega1
  gamma[21,16,] <- s * fid[8,3] * omega1
  gamma[21,17,] <- s * fid[8,4] * omega1
  gamma[21,18,] <- s * fid[8,5] * omega1
  gamma[21,19,] <- s * fid[8,6] * omega1
  gamma[21,20,] <- s * fid[8,7] * omega1
  gamma[21,21,] <- s * fid[8,8] * omega1
  gamma[21,22,] <- s * fid[8,9] * omega1
  gamma[21,23,] <- s * fid[8,10] * omega1
  gamma[21,24,] <- s * fid[8,11] * omega1
  gamma[21,25,] <- s * fid[8,12] * omega1
  gamma[21,26,] <- 1-s
  
  gamma[22,2,] <- s * fid[9,1] * (1-omega1)
  gamma[22,3,] <- s * fid[9,2] * (1-omega1)
  gamma[22,4,] <- s * fid[9,3] * (1-omega1)
  gamma[22,5,] <- s * fid[9,4] * (1-omega1)
  gamma[22,6,] <- s * fid[9,5] * (1-omega1)
  gamma[22,7,] <- s * fid[9,6] * (1-omega1)
  gamma[22,8,] <- s * fid[9,7] * (1-omega1)
  gamma[22,9,] <- s * fid[9,8] * (1-omega1)
  gamma[22,10,] <- s * fid[9,9] * (1-omega1)
  gamma[22,11,] <- s * fid[9,10] * (1-omega1)
  gamma[22,12,] <- s * fid[9,11] * (1-omega1)
  gamma[22,13,] <- s * fid[9,12] * (1-omega1)
  gamma[22,14,] <- s * fid[9,1] * omega1
  gamma[22,15,] <- s * fid[9,2] * omega1
  gamma[22,16,] <- s * fid[9,3] * omega1
  gamma[22,17,] <- s * fid[9,4] * omega1
  gamma[22,18,] <- s * fid[9,5] * omega1
  gamma[22,19,] <- s * fid[9,6] * omega1
  gamma[22,20,] <- s * fid[9,7] * omega1
  gamma[22,21,] <- s * fid[9,8] * omega1
  gamma[22,22,] <- s * fid[9,9] * omega1
  gamma[22,23,] <- s * fid[9,10] * omega1
  gamma[22,24,] <- s * fid[9,11] * omega1
  gamma[22,25,] <- s * fid[9,12] * omega1
  gamma[22,26,] <- 1-s
  
  gamma[23,2,] <- s * fid[10,1] * (1-omega1)
  gamma[23,3,] <- s * fid[10,2] * (1-omega1)
  gamma[23,4,] <- s * fid[10,3] * (1-omega1)
  gamma[23,5,] <- s * fid[10,4] * (1-omega1)
  gamma[23,6,] <- s * fid[10,5] * (1-omega1)
  gamma[23,7,] <- s * fid[10,6] * (1-omega1)
  gamma[23,8,] <- s * fid[10,7] * (1-omega1)
  gamma[23,9,] <- s * fid[10,8] * (1-omega1)
  gamma[23,10,] <- s * fid[10,9] * (1-omega1)
  gamma[23,11,] <- s * fid[10,10] * (1-omega1)
  gamma[23,12,] <- s * fid[10,11] * (1-omega1)
  gamma[23,13,] <- s * fid[10,12] * (1-omega1)
  gamma[23,14,] <- s * fid[10,1] * omega1
  gamma[23,15,] <- s * fid[10,2] * omega1
  gamma[23,16,] <- s * fid[10,3] * omega1
  gamma[23,17,] <- s * fid[10,4] * omega1
  gamma[23,18,] <- s * fid[10,5] * omega1
  gamma[23,19,] <- s * fid[10,6] * omega1
  gamma[23,20,] <- s * fid[10,7] * omega1
  gamma[23,21,] <- s * fid[10,8] * omega1
  gamma[23,22,] <- s * fid[10,9] * omega1
  gamma[23,23,] <- s * fid[10,10] * omega1
  gamma[23,24,] <- s * fid[10,11] * omega1
  gamma[23,25,] <- s * fid[10,12] * omega1
  gamma[23,26,] <- 1-s
  
  gamma[24,2,] <- s * fid[11,1] * (1-omega1)
  gamma[24,3,] <- s * fid[11,2] * (1-omega1)
  gamma[24,4,] <- s * fid[11,3] * (1-omega1)
  gamma[24,5,] <- s * fid[11,4] * (1-omega1)
  gamma[24,6,] <- s * fid[11,5] * (1-omega1)
  gamma[24,7,] <- s * fid[11,6] * (1-omega1)
  gamma[24,8,] <- s * fid[11,7] * (1-omega1)
  gamma[24,9,] <- s * fid[11,8] * (1-omega1)
  gamma[24,10,] <- s * fid[11,9] * (1-omega1)
  gamma[24,11,] <- s * fid[11,10] * (1-omega1)
  gamma[24,12,] <- s * fid[11,11] * (1-omega1)
  gamma[24,13,] <- s * fid[11,12] * (1-omega1)
  gamma[24,14,] <- s * fid[11,1] * omega1
  gamma[24,15,] <- s * fid[11,2] * omega1
  gamma[24,16,] <- s * fid[11,3] * omega1
  gamma[24,17,] <- s * fid[11,4] * omega1
  gamma[24,18,] <- s * fid[11,5] * omega1
  gamma[24,19,] <- s * fid[11,6] * omega1
  gamma[24,20,] <- s * fid[11,7] * omega1
  gamma[24,21,] <- s * fid[11,8] * omega1
  gamma[24,22,] <- s * fid[11,9] * omega1
  gamma[24,23,] <- s * fid[11,10] * omega1
  gamma[24,24,] <- s * fid[11,11] * omega1
  gamma[24,25,] <- s * fid[11,12] * omega1
  gamma[24,26,] <- 1-s
  
  gamma[25,2,] <- s * fid[12,1] * (1-omega1)
  gamma[25,3,] <- s * fid[12,2] * (1-omega1)
  gamma[25,4,] <- s * fid[12,3] * (1-omega1)
  gamma[25,5,] <- s * fid[12,4] * (1-omega1)
  gamma[25,6,] <- s * fid[12,5] * (1-omega1)
  gamma[25,7,] <- s * fid[12,6] * (1-omega1)
  gamma[25,8,] <- s * fid[12,7] * (1-omega1)
  gamma[25,9,] <- s * fid[12,8] * (1-omega1)
  gamma[25,10,] <- s * fid[12,9] * (1-omega1)
  gamma[25,11,] <- s * fid[12,10] * (1-omega1)
  gamma[25,12,] <- s * fid[12,11] * (1-omega1)
  gamma[25,13,] <- s * fid[12,12] * (1-omega1)
  gamma[25,14,] <- s * fid[12,1] * omega1
  gamma[25,15,] <- s * fid[12,2] * omega1
  gamma[25,16,] <- s * fid[12,3] * omega1
  gamma[25,17,] <- s * fid[12,4] * omega1
  gamma[25,18,] <- s * fid[12,5] * omega1
  gamma[25,19,] <- s * fid[12,6] * omega1
  gamma[25,20,] <- s * fid[12,7] * omega1
  gamma[25,21,] <- s * fid[12,8] * omega1
  gamma[25,22,] <- s * fid[12,9] * omega1
  gamma[25,23,] <- s * fid[12,10] * omega1
  gamma[25,24,] <- s * fid[12,11] * omega1
  gamma[25,25,] <- s * fid[12,12] * omega1
  gamma[25,26,] <- 1-s
  
  gamma[26,26,] <- 1
  
  # transition probability matrices within periods
  gammat <- NULL
  for (t in 1:T)  {
    gammat[[t]] <- array(0, dim=c(3, 3, K[t]-1))
    gammat[[t]][1,1,] <- 1 - betastar[[t]][2:K[t]]
    gammat[[t]][1,2,] <- betastar[[t]][2:K[t]]
    
    gammat[[t]][2,2,] <- phistar[[t]]
    gammat[[t]][2,3,] <- 1-phistar[[t]]
    
    gammat[[t]][3,3,] <- 1
  }
  
    
  omega <- matrix(0, 2, 1)
  omega[1,1] <- omega1
  omega[2,1] <- (1-omega1)

  
  # observation matrices within periods
  pmat <- NULL
  for (t in 1:T)  {
    pmat[[t]] <- array(0, dim=c(3,3,2,K[t]))
    pmat[[t]][1,1,1,] <- 1
    pmat[[t]][2,2,1,] <- 1-p[[t]]
    pmat[[t]][3,3,1,] <- 1
    
    pmat[[t]][2,2,2,] <- p[[t]]
  }
  
  # state probabilities for period 1
  pione <- NULL
  pione[1] <- 1-rstar[1]
  pione[2] <- rstar[1]*alpha[1]
  pione[3] <- rstar[1]*alpha[2]
  pione[4] <- rstar[1]*alpha[3]
  pione[5] <- rstar[1]*alpha[4]
  pione[6] <- rstar[1]*alpha[5]
  pione[7] <- rstar[1]*alpha[6]
  pione[8] <- rstar[1]*alpha[7]
  pione[9] <- rstar[1]*alpha[8]
  pione[10] <- rstar[1]*alpha[9]
  pione[11] <- rstar[1]*alpha[10]
  pione[12] <- rstar[1]*alpha[11]
  pione[13] <- rstar[1]*alpha[12]
  pione[14] <- 0
  pione[15] <- 0
  pione[16] <- 0
  pione[17] <- 0
  pione[18] <- 0
  pione[19] <- 0
  pione[20] <- 0
  pione[21] <- 0
  pione[22] <- 0
  pione[23] <- 0
  pione[24] <- 0
  pione[25] <- 0
  pione[26] <- 0
  
  # state probabilities for occasion 1 within periods
  pionet <- NULL
  for (t in 1:T)  {
    pionet[[t]] <- rep(0,3)
    pionet[[t]][1] <- 1-betastar[[t]][1]
    pionet[[t]][2] <- betastar[[t]][1]
  }
  
  # return all parameters
  return(list('n.missed'=n.missed, 'N'=N, 'r'=r, 'rstar'=rstar, 's'=s, 'betaintercept'=betaintercept, 'betagradient'=betagradient, 'betalogit'=betalogit, 'betastar'=betastar, 'omega1'=omega1, 'omega2' = omega2, 'p'=p, 'sigma.primary' = sigma.primary, 'fid' = fid, 'phiintercept'=phiintercept, 'phigradient'=phigradient, 'philogit'=philogit, 'phistar'=phistar, 'gamma'=gamma, 'gammat'=gammat, 'pmat'=pmat, 'pione'=pione, 'pionet'=pionet))
}

# -------------------------------------------------------------------------------
# HMM likelihood
# -------------------------------------------------------------------------------

# Name: likelihood
# Objective: To calculate the negative log-likelihood of the HMM model for the greated crested newt data
# Inputs: param - vector of parameters to be passed to param.unpack
#         n - number of observed individuals
#         T - number of years
#         K - number of occasions in each year
#         X - list of capture histories
#         Z - attendance history
#         n.hist - frequency of histories
# Outputs: lik - value of the negative log-likelihood

likelihood <- function(param, n, T, K, X, Z, n.hist, distSq)  {
  
  # occasions that mark the start of each period
  occ <- rep(1, T+1)
  for (t in 2:(T+1))  {
    occ[t] <- sum(K[1:(t-1)])
  }
  
  # unpack parameter vector
  params <- param.unpack(param, n, T, K, distSq)
  
  # multinomial term
  term1 <- (lfactorial(params$N) - lfactorial(params$n.missed) - sum(lfactorial(n.hist)))
  
  # unobserved individuals within periods
  likzerot <- rep(0, T)
  for (t in 1:T)  {
    test <- params$pionet[[t]]
    test <- test%*%params$pmat[[t]][,,1,1]
    for (k in 1:(K[t]-1))  {
      test <- test%*%params$gammat[[t]][,,k]
      test <- test%*%params$pmat[[t]][,,1,k+1]
    }
    test <- test%*%matrix(1, nrow=3, ncol=1)
    likzerot[t] <- test
  }
  
  # observation matrix for missed individuals
  pzero <- array(0, dim=c(26,26,T))
  pzero[1,1,] <- 1
  pzero[2,2,] <- likzerot
  pzero[3,3,] <- likzerot
  pzero[4,4,] <- likzerot
  pzero[5,5,] <- likzerot
  pzero[6,6,] <- likzerot
  pzero[7,7,] <- likzerot
  pzero[8,8,] <- likzerot
  pzero[9,9,] <- likzerot
  pzero[10,10,] <- likzerot
  pzero[11,11,] <- likzerot
  pzero[12,12,] <- likzerot
  pzero[13,13,] <- likzerot
  pzero[14,14,] <- 1
  pzero[15,15,] <- 1
  pzero[16,16,] <- 1
  pzero[17,17,] <- 1
  pzero[18,18,] <- 1
  pzero[19,19,] <- 1
  pzero[20,20,] <- 1
  pzero[21,21,] <- 1
  pzero[22,22,] <- 1
  pzero[23,23,] <- 1
  pzero[24,24,] <- 1
  pzero[25,25,] <- 1
  pzero[26,26,] <- 1
  
  # unobserved individuals 
  likzero <- params$pione
  likzero <- likzero%*%pzero[,,1]
  for (t in 1:(T-1))  {
    likzero <- likzero%*%params$gamma[,,t]
    likzero <- likzero%*%pzero[,,t+1]
  }
  likzero <- likzero%*%matrix(1, nrow=26, ncol=1)
  term2 <- (params$n.missed)*log(likzero)
  
  # observed individuals
  likobs <- rep(0,n)
  for (i in 1:n)  {
    likobst <- rep(0,T)
    for (t in 1:T)  {
      if (Z[i,t] == 0)  {
        likobst[t] <- likzerot[t]
      } else if (Z[i,t] == 1 | Z[i,t] == 2 | Z[i,t] == 3 | Z[i,t] == 4 | Z[i,t] == 5 | Z[i,t] == 6 |
                 Z[i,t] == 7 | Z[i,t] == 8 | Z[i,t] == 9 | Z[i,t] == 10 |Z[i,t] == 11 | Z[i,t] == 12)  {
        test <- params$pionet[[t]]
        test <- test%*%params$pmat[[t]][,,X[[t]][i,1]+1,1]
        for (k in 1:(K[t]-1))  {
          test <- test%*%params$gammat[[t]][,,k]
          test <- test%*%params$pmat[[t]][,,X[[t]][i,k+1]+1,k+1]
          }
        likobst[t] <- test%*%matrix(1, nrow=3, ncol=1)
      }
    }
    
    # observation matrices
    pobs <- array(0, dim=c(26,26,13,T))
    pobs[1,1,1,] <- 1
    pobs[2,2,1,] <- likzerot
    pobs[3,3,1,] <- likzerot
    pobs[4,4,1,] <- likzerot
    pobs[5,5,1,] <- likzerot
    pobs[6,6,1,] <- likzerot
    pobs[7,7,1,] <- likzerot
    pobs[8,8,1,] <- likzerot
    pobs[9,9,1,] <- likzerot
    pobs[10,10,1,] <- likzerot
    pobs[11,11,1,] <- likzerot
    pobs[12,12,1,] <- likzerot
    pobs[13,13,1,] <- likzerot
    pobs[14,14,1,] <- 1
    pobs[14,14,1,] <- 1
    pobs[15,15,1,] <- 1
    pobs[16,16,1,] <- 1
    pobs[17,17,1,] <- 1
    pobs[18,18,1,] <- 1
    pobs[19,19,1,] <- 1
    pobs[20,20,1,] <- 1
    pobs[21,21,1,] <- 1
    pobs[22,22,1,] <- 1
    pobs[23,23,1,] <- 1
    pobs[24,24,1,] <- 1
    pobs[25,25,1,] <- 1
    pobs[26,26,1,] <- 1
    
    pobs[2,2,2,] <- likobst
    pobs[3,3,3,] <- likobst
    pobs[4,4,4,] <- likobst
    pobs[5,5,5,] <- likobst
    pobs[6,6,6,] <- likobst
    pobs[7,7,7,] <- likobst
    pobs[8,8,8,] <- likobst
    pobs[9,9,9,] <- likobst
    pobs[10,10,10,] <- likobst
    pobs[11,11,11,] <- likobst
    pobs[12,12,12,] <- likobst
    pobs[13,13,13,] <- likobst
    
    
    # likelihood contribution
    test <- params$pione
    test <- test%*%pobs[,,Z[i,1]+1,1]
    for (t in 1:(T-1))  {
      test <- test%*%params$gamma[,,t]
      test <- test%*%pobs[,,Z[i,t+1]+1,t+1]
    }
    likobs[i] <- test%*%matrix(1, nrow=26, ncol=1)
  }
  
  term3 <- sum(log(likobs))
  
  # full log-likelihood
  lik <- term1+term2+term3
  
  # negate likelihood
  lik <- -lik
  
  if (isTRUE(lik < 0))  {
    lik <- 10000
  }
  
  # return the value of the likelihood
  #print(lik)
  #print(param)
  return(lik)
}
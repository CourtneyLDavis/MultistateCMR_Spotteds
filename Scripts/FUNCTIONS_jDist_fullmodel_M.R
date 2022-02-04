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
# This script defines the function to unpack the parameter vector for the male
# model that estimates inter- and intra-annual movement as a function of Euclidean 
# distance between sites.
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
#               - [32] - sigma - scale parameter for rate of decay of within-period dispersal
#               - [33:44] - within-period site-specific fidelity (probability of no movement)
#               - [45] - intercept for logistic regression on retention (constant) 
#               - [46:51] - log gradient for logistic regression on retention
#               - [52] - sigma.primary - scale parameter for rate of decay of between-period dispersal
#               - [53:64] - between-period site-specific fidelity (probability of no movement)
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
#          alpha - initial discrete state probabilities
#          p - capture probabilities
#          sigma - scale parameter that determines the rate of decay within-period dispersal probability as a function of distance
#          sigma.primary - scale parameter that determines the rate of decay between-period dispersal probability as a function of distance
#          psi - movement probabilities within periods
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
  
  # observable states transition probabilities
  sigma <- exp(param[32])

  # movement as a function of Euclidean distance between sites
  move <- array(0, dim = c(12,12))
  move[1,2] <- exp((-distSq[1,2]^2)/(2*sigma^2))
  move[1,3] <- exp((-distSq[1,3]^2)/(2*sigma^2))
  move[1,4] <- exp((-distSq[1,4]^2)/(2*sigma^2))
  move[1,5] <- exp((-distSq[1,5]^2)/(2*sigma^2))
  move[1,6] <- exp((-distSq[1,6]^2)/(2*sigma^2))
  move[1,7] <- exp((-distSq[1,7]^2)/(2*sigma^2))
  move[1,8] <- exp((-distSq[1,8]^2)/(2*sigma^2))
  move[1,9] <- exp((-distSq[1,9]^2)/(2*sigma^2))
  move[1,10] <- exp((-distSq[1,10]^2)/(2*sigma^2))
  move[1,11] <- exp((-distSq[1,11]^2)/(2*sigma^2))
  move[1,12] <- exp((-distSq[1,12]^2)/(2*sigma^2))
  
  move[2,1] <- exp((-distSq[2,1]^2)/(2*sigma^2))
  move[2,3] <- exp((-distSq[2,3]^2)/(2*sigma^2))
  move[2,4] <- exp((-distSq[2,4]^2)/(2*sigma^2))
  move[2,5] <- exp((-distSq[2,5]^2)/(2*sigma^2))
  move[2,6] <- exp((-distSq[2,6]^2)/(2*sigma^2))
  move[2,7] <- exp((-distSq[2,7]^2)/(2*sigma^2))
  move[2,8] <- exp((-distSq[2,8]^2)/(2*sigma^2))
  move[2,9] <- exp((-distSq[2,9]^2)/(2*sigma^2))
  move[2,10] <- exp((-distSq[2,10]^2)/(2*sigma^2))
  move[2,11] <- exp((-distSq[2,11]^2)/(2*sigma^2))
  move[2,12] <- exp((-distSq[2,12]^2)/(2*sigma^2))
  
  move[3,1] <- exp((-distSq[3,1]^2)/(2*sigma^2))
  move[3,2] <- exp((-distSq[3,2]^2)/(2*sigma^2))
  move[3,4] <- exp((-distSq[3,4]^2)/(2*sigma^2))
  move[3,5] <- exp((-distSq[3,5]^2)/(2*sigma^2))
  move[3,6] <- exp((-distSq[3,6]^2)/(2*sigma^2))
  move[3,7] <- exp((-distSq[3,7]^2)/(2*sigma^2))
  move[3,8] <- exp((-distSq[3,8]^2)/(2*sigma^2))
  move[3,9] <- exp((-distSq[3,9]^2)/(2*sigma^2))
  move[3,10] <- exp((-distSq[3,10]^2)/(2*sigma^2))
  move[3,11] <- exp((-distSq[3,11]^2)/(2*sigma^2))
  move[3,12] <- exp((-distSq[3,12]^2)/(2*sigma^2))
  
  move[4,1] <- exp((-distSq[4,1]^2)/(2*sigma^2))
  move[4,2] <- exp((-distSq[4,2]^2)/(2*sigma^2))
  move[4,3] <- exp((-distSq[4,3]^2)/(2*sigma^2))
  move[4,5] <- exp((-distSq[4,5]^2)/(2*sigma^2))
  move[4,6] <- exp((-distSq[4,6]^2)/(2*sigma^2))
  move[4,7] <- exp((-distSq[4,7]^2)/(2*sigma^2))
  move[4,8] <- exp((-distSq[4,8]^2)/(2*sigma^2))
  move[4,9] <- exp((-distSq[4,9]^2)/(2*sigma^2))
  move[4,10] <- exp((-distSq[4,10]^2)/(2*sigma^2))
  move[4,11] <- exp((-distSq[4,11]^2)/(2*sigma^2))
  move[4,12] <- exp((-distSq[4,12]^2)/(2*sigma^2))
  
  move[5,1] <- exp((-distSq[5,1]^2)/(2*sigma^2))
  move[5,2] <- exp((-distSq[5,2]^2)/(2*sigma^2))
  move[5,3] <- exp((-distSq[5,3]^2)/(2*sigma^2))
  move[5,4] <- exp((-distSq[5,4]^2)/(2*sigma^2))
  move[5,6] <- exp((-distSq[5,6]^2)/(2*sigma^2))
  move[5,7] <- exp((-distSq[5,7]^2)/(2*sigma^2))
  move[5,8] <- exp((-distSq[5,8]^2)/(2*sigma^2))
  move[5,9] <- exp((-distSq[5,9]^2)/(2*sigma^2))
  move[5,10] <- exp((-distSq[5,10]^2)/(2*sigma^2))
  move[5,11] <- exp((-distSq[5,11]^2)/(2*sigma^2))
  move[5,12] <- exp((-distSq[5,12]^2)/(2*sigma^2))
  
  move[6,1] <- exp((-distSq[6,1]^2)/(2*sigma^2))
  move[6,2] <- exp((-distSq[6,2]^2)/(2*sigma^2))
  move[6,3] <- exp((-distSq[6,3]^2)/(2*sigma^2))
  move[6,4] <- exp((-distSq[6,4]^2)/(2*sigma^2))
  move[6,5] <- exp((-distSq[6,5]^2)/(2*sigma^2))
  move[6,7] <- exp((-distSq[6,7]^2)/(2*sigma^2))
  move[6,8] <- exp((-distSq[6,8]^2)/(2*sigma^2))
  move[6,9] <- exp((-distSq[6,9]^2)/(2*sigma^2))
  move[6,10] <- exp((-distSq[6,10]^2)/(2*sigma^2))
  move[6,11] <- exp((-distSq[6,11]^2)/(2*sigma^2))
  move[6,12] <- exp((-distSq[6,12]^2)/(2*sigma^2))
  
  move[7,1] <- exp((-distSq[7,1]^2)/(2*sigma^2))
  move[7,2] <- exp((-distSq[7,2]^2)/(2*sigma^2))
  move[7,3] <- exp((-distSq[7,3]^2)/(2*sigma^2))
  move[7,4] <- exp((-distSq[7,4]^2)/(2*sigma^2))
  move[7,5] <- exp((-distSq[7,5]^2)/(2*sigma^2))
  move[7,6] <- exp((-distSq[7,6]^2)/(2*sigma^2))
  move[7,8] <- exp((-distSq[7,8]^2)/(2*sigma^2))
  move[7,9] <- exp((-distSq[7,9]^2)/(2*sigma^2))
  move[7,10] <- exp((-distSq[7,10]^2)/(2*sigma^2))
  move[7,11] <- exp((-distSq[7,11]^2)/(2*sigma^2))
  move[7,12] <- exp((-distSq[7,12]^2)/(2*sigma^2))
  
  move[8,1] <- exp((-distSq[8,1]^2)/(2*sigma^2))
  move[8,2] <- exp((-distSq[8,2]^2)/(2*sigma^2))
  move[8,3] <- exp((-distSq[8,3]^2)/(2*sigma^2))
  move[8,4] <- exp((-distSq[8,4]^2)/(2*sigma^2))
  move[8,5] <- exp((-distSq[8,5]^2)/(2*sigma^2))
  move[8,6] <- exp((-distSq[8,6]^2)/(2*sigma^2))
  move[8,7] <- exp((-distSq[8,7]^2)/(2*sigma^2))
  move[8,9] <- exp((-distSq[8,9]^2)/(2*sigma^2))
  move[8,10] <- exp((-distSq[8,10]^2)/(2*sigma^2))
  move[8,11] <- exp((-distSq[8,11]^2)/(2*sigma^2))
  move[8,12] <- exp((-distSq[8,12]^2)/(2*sigma^2))
  
  move[9,1] <- exp((-distSq[9,1]^2)/(2*sigma^2))
  move[9,2] <- exp((-distSq[9,2]^2)/(2*sigma^2))
  move[9,3] <- exp((-distSq[9,3]^2)/(2*sigma^2))
  move[9,4] <- exp((-distSq[9,4]^2)/(2*sigma^2))
  move[9,5] <- exp((-distSq[9,5]^2)/(2*sigma^2))
  move[9,6] <- exp((-distSq[9,6]^2)/(2*sigma^2))
  move[9,7] <- exp((-distSq[9,7]^2)/(2*sigma^2))
  move[9,8] <- exp((-distSq[9,8]^2)/(2*sigma^2))
  move[9,10] <- exp((-distSq[9,10]^2)/(2*sigma^2))
  move[9,11] <- exp((-distSq[9,11]^2)/(2*sigma^2))
  move[9,12] <- exp((-distSq[9,12]^2)/(2*sigma^2))
  
  move[10,1] <- exp((-distSq[10,1]^2)/(2*sigma^2))
  move[10,2] <- exp((-distSq[10,2]^2)/(2*sigma^2))
  move[10,3] <- exp((-distSq[10,3]^2)/(2*sigma^2))
  move[10,4] <- exp((-distSq[10,4]^2)/(2*sigma^2))
  move[10,5] <- exp((-distSq[10,5]^2)/(2*sigma^2))
  move[10,6] <- exp((-distSq[10,6]^2)/(2*sigma^2))
  move[10,7] <- exp((-distSq[10,7]^2)/(2*sigma^2))
  move[10,8] <- exp((-distSq[10,8]^2)/(2*sigma^2))
  move[10,9] <- exp((-distSq[10,9]^2)/(2*sigma^2))
  move[10,11] <- exp((-distSq[10,11]^2)/(2*sigma^2))
  move[10,12] <- exp((-distSq[10,12]^2)/(2*sigma^2))
  
  move[11,1] <- exp((-distSq[11,1]^2)/(2*sigma^2))
  move[11,2] <- exp((-distSq[11,2]^2)/(2*sigma^2))
  move[11,3] <- exp((-distSq[11,3]^2)/(2*sigma^2))
  move[11,4] <- exp((-distSq[11,4]^2)/(2*sigma^2))
  move[11,5] <- exp((-distSq[11,5]^2)/(2*sigma^2))
  move[11,6] <- exp((-distSq[11,6]^2)/(2*sigma^2))
  move[11,7] <- exp((-distSq[11,7]^2)/(2*sigma^2))
  move[11,8] <- exp((-distSq[11,8]^2)/(2*sigma^2))
  move[11,9] <- exp((-distSq[11,9]^2)/(2*sigma^2))
  move[11,10] <- exp((-distSq[11,10]^2)/(2*sigma^2))
  move[11,12] <- exp((-distSq[11,12]^2)/(2*sigma^2))
  
  move[12,1] <- exp((-distSq[12,1]^2)/(2*sigma^2))
  move[12,2] <- exp((-distSq[12,2]^2)/(2*sigma^2))
  move[12,3] <- exp((-distSq[12,3]^2)/(2*sigma^2))
  move[12,4] <- exp((-distSq[12,4]^2)/(2*sigma^2))
  move[12,5] <- exp((-distSq[12,5]^2)/(2*sigma^2))
  move[12,6] <- exp((-distSq[12,6]^2)/(2*sigma^2))
  move[12,7] <- exp((-distSq[12,7]^2)/(2*sigma^2))
  move[12,8] <- exp((-distSq[12,8]^2)/(2*sigma^2))
  move[12,9] <- exp((-distSq[12,9]^2)/(2*sigma^2))
  move[12,10] <- exp((-distSq[12,10]^2)/(2*sigma^2))
  move[12,11] <- exp((-distSq[12,11]^2)/(2*sigma^2))
  
  beta.stay.secondary <- matrix(0, 12, 1) # probability an animal stays in a pond within a year
  beta.stay.secondary <- 1/(1+exp(-(param[33:44])))

  psi <- array(0, dim=c(12,12))
  psi[1,2] <- (1-beta.stay.secondary[1])*(move[1,2]/(move[1,2]+move[1,3]+move[1,4]+move[1,5]+move[1,6]+move[1,7]+move[1,8]+move[1,9]+move[1,10]+move[1,11]+move[1,12]))
  psi[1,3] <- (1-beta.stay.secondary[1])*(move[1,3]/(move[1,2]+move[1,3]+move[1,4]+move[1,5]+move[1,6]+move[1,7]+move[1,8]+move[1,9]+move[1,10]+move[1,11]+move[1,12]))
  psi[1,4] <- (1-beta.stay.secondary[1])*(move[1,4]/(move[1,2]+move[1,3]+move[1,4]+move[1,5]+move[1,6]+move[1,7]+move[1,8]+move[1,9]+move[1,10]+move[1,11]+move[1,12]))
  psi[1,5] <- (1-beta.stay.secondary[1])*(move[1,5]/(move[1,2]+move[1,3]+move[1,4]+move[1,5]+move[1,6]+move[1,7]+move[1,8]+move[1,9]+move[1,10]+move[1,11]+move[1,12]))
  psi[1,6] <- (1-beta.stay.secondary[1])*(move[1,6]/(move[1,2]+move[1,3]+move[1,4]+move[1,5]+move[1,6]+move[1,7]+move[1,8]+move[1,9]+move[1,10]+move[1,11]+move[1,12]))
  psi[1,7] <- (1-beta.stay.secondary[1])*(move[1,7]/(move[1,2]+move[1,3]+move[1,4]+move[1,5]+move[1,6]+move[1,7]+move[1,8]+move[1,9]+move[1,10]+move[1,11]+move[1,12]))
  psi[1,8] <- (1-beta.stay.secondary[1])*(move[1,8]/(move[1,2]+move[1,3]+move[1,4]+move[1,5]+move[1,6]+move[1,7]+move[1,8]+move[1,9]+move[1,10]+move[1,11]+move[1,12]))
  psi[1,9] <- (1-beta.stay.secondary[1])*(move[1,9]/(move[1,2]+move[1,3]+move[1,4]+move[1,5]+move[1,6]+move[1,7]+move[1,8]+move[1,9]+move[1,10]+move[1,11]+move[1,12]))
  psi[1,10] <- (1-beta.stay.secondary[1])*(move[1,10]/(move[1,2]+move[1,3]+move[1,4]+move[1,5]+move[1,6]+move[1,7]+move[1,8]+move[1,9]+move[1,10]+move[1,11]+move[1,12]))
  psi[1,11] <- (1-beta.stay.secondary[1])*(move[1,11]/(move[1,2]+move[1,3]+move[1,4]+move[1,5]+move[1,6]+move[1,7]+move[1,8]+move[1,9]+move[1,10]+move[1,11]+move[1,12]))
  psi[1,12] <- (1-beta.stay.secondary[1])*(move[1,12]/(move[1,2]+move[1,3]+move[1,4]+move[1,5]+move[1,6]+move[1,7]+move[1,8]+move[1,9]+move[1,10]+move[1,11]+move[1,12]))
  psi[1,1] <- beta.stay.secondary[1]
  
  psi[2,1] <- (1-beta.stay.secondary[2])*(move[2,1]/(move[2,1]+move[2,3]+move[2,4]+move[2,5]+move[2,6]+move[2,7]+move[2,8]+move[2,9]+move[2,10]+move[2,11]+move[2,12]))
  psi[2,3] <- (1-beta.stay.secondary[2])*(move[2,3]/(move[2,1]+move[2,3]+move[2,4]+move[2,5]+move[2,6]+move[2,7]+move[2,8]+move[2,9]+move[2,10]+move[2,11]+move[2,12]))
  psi[2,4] <- (1-beta.stay.secondary[2])*(move[2,4]/(move[2,1]+move[2,3]+move[2,4]+move[2,5]+move[2,6]+move[2,7]+move[2,8]+move[2,9]+move[2,10]+move[2,11]+move[2,12]))
  psi[2,5] <- (1-beta.stay.secondary[2])*(move[2,5]/(move[2,1]+move[2,3]+move[2,4]+move[2,5]+move[2,6]+move[2,7]+move[2,8]+move[2,9]+move[2,10]+move[2,11]+move[2,12]))
  psi[2,6] <- (1-beta.stay.secondary[2])*(move[2,6]/(move[2,1]+move[2,3]+move[2,4]+move[2,5]+move[2,6]+move[2,7]+move[2,8]+move[2,9]+move[2,10]+move[2,11]+move[2,12]))
  psi[2,7] <- (1-beta.stay.secondary[2])*(move[2,7]/(move[2,1]+move[2,3]+move[2,4]+move[2,5]+move[2,6]+move[2,7]+move[2,8]+move[2,9]+move[2,10]+move[2,11]+move[2,12]))
  psi[2,8] <- (1-beta.stay.secondary[2])*(move[2,8]/(move[2,1]+move[2,3]+move[2,4]+move[2,5]+move[2,6]+move[2,7]+move[2,8]+move[2,9]+move[2,10]+move[2,11]+move[2,12]))
  psi[2,9] <- (1-beta.stay.secondary[2])*(move[2,9]/(move[2,1]+move[2,3]+move[2,4]+move[2,5]+move[2,6]+move[2,7]+move[2,8]+move[2,9]+move[2,10]+move[2,11]+move[2,12]))
  psi[2,10] <- (1-beta.stay.secondary[2])*(move[2,10]/(move[2,1]+move[2,3]+move[2,4]+move[2,5]+move[2,6]+move[2,7]+move[2,8]+move[2,9]+move[2,10]+move[2,11]+move[2,12]))
  psi[2,11] <- (1-beta.stay.secondary[2])*(move[2,11]/(move[2,1]+move[2,3]+move[2,4]+move[2,5]+move[2,6]+move[2,7]+move[2,8]+move[2,9]+move[2,10]+move[2,11]+move[2,12]))
  psi[2,12] <- (1-beta.stay.secondary[2])*(move[2,12]/(move[2,1]+move[2,3]+move[2,4]+move[2,5]+move[2,6]+move[2,7]+move[2,8]+move[2,9]+move[2,10]+move[2,11]+move[2,12]))
  psi[2,2] <- beta.stay.secondary[2]
  
  psi[3,1] <- (1-beta.stay.secondary[3])*(move[3,1]/(move[3,1]+move[3,2]+move[3,4]+move[3,5]+move[3,6]+move[3,7]+move[3,8]+move[3,9]+move[3,10]+move[3,11]+move[3,12]))
  psi[3,2] <- (1-beta.stay.secondary[3])*(move[3,2]/(move[3,1]+move[3,2]+move[3,4]+move[3,5]+move[3,6]+move[3,7]+move[3,8]+move[3,9]+move[3,10]+move[3,11]+move[3,12]))
  psi[3,4] <- (1-beta.stay.secondary[3])*(move[3,4]/(move[3,1]+move[3,2]+move[3,4]+move[3,5]+move[3,6]+move[3,7]+move[3,8]+move[3,9]+move[3,10]+move[3,11]+move[3,12]))
  psi[3,5] <- (1-beta.stay.secondary[3])*(move[3,5]/(move[3,1]+move[3,2]+move[3,4]+move[3,5]+move[3,6]+move[3,7]+move[3,8]+move[3,9]+move[3,10]+move[3,11]+move[3,12]))
  psi[3,6] <- (1-beta.stay.secondary[3])*(move[3,6]/(move[3,1]+move[3,2]+move[3,4]+move[3,5]+move[3,6]+move[3,7]+move[3,8]+move[3,9]+move[3,10]+move[3,11]+move[3,12]))
  psi[3,7] <- (1-beta.stay.secondary[3])*(move[3,7]/(move[3,1]+move[3,2]+move[3,4]+move[3,5]+move[3,6]+move[3,7]+move[3,8]+move[3,9]+move[3,10]+move[3,11]+move[3,12]))
  psi[3,8] <- (1-beta.stay.secondary[3])*(move[3,8]/(move[3,1]+move[3,2]+move[3,4]+move[3,5]+move[3,6]+move[3,7]+move[3,8]+move[3,9]+move[3,10]+move[3,11]+move[3,12]))
  psi[3,9] <- (1-beta.stay.secondary[3])*(move[3,9]/(move[3,1]+move[3,2]+move[3,4]+move[3,5]+move[3,6]+move[3,7]+move[3,8]+move[3,9]+move[3,10]+move[3,11]+move[3,12]))
  psi[3,10] <- (1-beta.stay.secondary[3])*(move[3,10]/(move[3,1]+move[3,2]+move[3,4]+move[3,5]+move[3,6]+move[3,7]+move[3,8]+move[3,9]+move[3,10]+move[3,11]+move[3,12]))
  psi[3,11] <- (1-beta.stay.secondary[3])*(move[3,11]/(move[3,1]+move[3,2]+move[3,4]+move[3,5]+move[3,6]+move[3,7]+move[3,8]+move[3,9]+move[3,10]+move[3,11]+move[3,12]))
  psi[3,12] <- (1-beta.stay.secondary[3])*(move[3,12]/(move[3,1]+move[3,2]+move[3,4]+move[3,5]+move[3,6]+move[3,7]+move[3,8]+move[3,9]+move[3,10]+move[3,11]+move[3,12]))
  psi[3,3] <- beta.stay.secondary[3]
  
  psi[4,1] <- (1-beta.stay.secondary[4])*(move[4,1]/(move[4,1]+move[4,2]+move[4,3]+move[4,5]+move[4,6]+move[4,7]+move[4,8]+move[4,9]+move[4,10]+move[4,11]+move[4,12]))
  psi[4,2] <- (1-beta.stay.secondary[4])*(move[4,2]/(move[4,1]+move[4,2]+move[4,3]+move[4,5]+move[4,6]+move[4,7]+move[4,8]+move[4,9]+move[4,10]+move[4,11]+move[4,12]))
  psi[4,3] <- (1-beta.stay.secondary[4])*(move[4,3]/(move[4,1]+move[4,2]+move[4,3]+move[4,5]+move[4,6]+move[4,7]+move[4,8]+move[4,9]+move[4,10]+move[4,11]+move[4,12]))
  psi[4,5] <- (1-beta.stay.secondary[4])*(move[4,5]/(move[4,1]+move[4,2]+move[4,3]+move[4,5]+move[4,6]+move[4,7]+move[4,8]+move[4,9]+move[4,10]+move[4,11]+move[4,12]))
  psi[4,6] <- (1-beta.stay.secondary[4])*(move[4,6]/(move[4,1]+move[4,2]+move[4,3]+move[4,5]+move[4,6]+move[4,7]+move[4,8]+move[4,9]+move[4,10]+move[4,11]+move[4,12]))
  psi[4,7] <- (1-beta.stay.secondary[4])*(move[4,7]/(move[4,1]+move[4,2]+move[4,3]+move[4,5]+move[4,6]+move[4,7]+move[4,8]+move[4,9]+move[4,10]+move[4,11]+move[4,12]))
  psi[4,8] <- (1-beta.stay.secondary[4])*(move[4,8]/(move[4,1]+move[4,2]+move[4,3]+move[4,5]+move[4,6]+move[4,7]+move[4,8]+move[4,9]+move[4,10]+move[4,11]+move[4,12]))
  psi[4,9] <- (1-beta.stay.secondary[4])*(move[4,9]/(move[4,1]+move[4,2]+move[4,3]+move[4,5]+move[4,6]+move[4,7]+move[4,8]+move[4,9]+move[4,10]+move[4,11]+move[4,12]))
  psi[4,10] <- (1-beta.stay.secondary[4])*(move[4,10]/(move[4,1]+move[4,2]+move[4,3]+move[4,5]+move[4,6]+move[4,7]+move[4,8]+move[4,9]+move[4,10]+move[4,11]+move[4,12]))
  psi[4,11] <- (1-beta.stay.secondary[4])*(move[4,11]/(move[4,1]+move[4,2]+move[4,3]+move[4,5]+move[4,6]+move[4,7]+move[4,8]+move[4,9]+move[4,10]+move[4,11]+move[4,12]))
  psi[4,12] <- (1-beta.stay.secondary[4])*(move[4,12]/(move[4,1]+move[4,2]+move[4,3]+move[4,5]+move[4,6]+move[4,7]+move[4,8]+move[4,9]+move[4,10]+move[4,11]+move[4,12]))
  psi[4,4] <- beta.stay.secondary[4]
  
  psi[5,1] <- (1-beta.stay.secondary[5])*(move[5,1]/(move[5,1]+move[5,2]+move[5,3]+move[5,4]+move[5,6]+move[5,7]+move[5,8]+move[5,9]+move[5,10]+move[5,11]+move[5,12]))
  psi[5,2] <- (1-beta.stay.secondary[5])*(move[5,2]/(move[5,1]+move[5,2]+move[5,3]+move[5,4]+move[5,6]+move[5,7]+move[5,8]+move[5,9]+move[5,10]+move[5,11]+move[5,12]))
  psi[5,3] <- (1-beta.stay.secondary[5])*(move[5,3]/(move[5,1]+move[5,2]+move[5,3]+move[5,4]+move[5,6]+move[5,7]+move[5,8]+move[5,9]+move[5,10]+move[5,11]+move[5,12]))
  psi[5,4] <- (1-beta.stay.secondary[5])*(move[5,4]/(move[5,1]+move[5,2]+move[5,3]+move[5,4]+move[5,6]+move[5,7]+move[5,8]+move[5,9]+move[5,10]+move[5,11]+move[5,12]))
  psi[5,6] <- (1-beta.stay.secondary[5])*(move[5,6]/(move[5,1]+move[5,2]+move[5,3]+move[5,4]+move[5,6]+move[5,7]+move[5,8]+move[5,9]+move[5,10]+move[5,11]+move[5,12]))
  psi[5,7] <- (1-beta.stay.secondary[5])*(move[5,7]/(move[5,1]+move[5,2]+move[5,3]+move[5,4]+move[5,6]+move[5,7]+move[5,8]+move[5,9]+move[5,10]+move[5,11]+move[5,12]))
  psi[5,8] <- (1-beta.stay.secondary[5])*(move[5,8]/(move[5,1]+move[5,2]+move[5,3]+move[5,4]+move[5,6]+move[5,7]+move[5,8]+move[5,9]+move[5,10]+move[5,11]+move[5,12]))
  psi[5,9] <- (1-beta.stay.secondary[5])*(move[5,9]/(move[5,1]+move[5,2]+move[5,3]+move[5,4]+move[5,6]+move[5,7]+move[5,8]+move[5,9]+move[5,10]+move[5,11]+move[5,12]))
  psi[5,10] <- (1-beta.stay.secondary[5])*(move[5,10]/(move[5,1]+move[5,2]+move[5,3]+move[5,4]+move[5,6]+move[5,7]+move[5,8]+move[5,9]+move[5,10]+move[5,11]+move[5,12]))
  psi[5,11] <- (1-beta.stay.secondary[5])*(move[5,11]/(move[5,1]+move[5,2]+move[5,3]+move[5,4]+move[5,6]+move[5,7]+move[5,8]+move[5,9]+move[5,10]+move[5,11]+move[5,12]))
  psi[5,12] <- (1-beta.stay.secondary[5])*(move[5,12]/(move[5,1]+move[5,2]+move[5,3]+move[5,4]+move[5,6]+move[5,7]+move[5,8]+move[5,9]+move[5,10]+move[5,11]+move[5,12]))
  psi[5,5] <- beta.stay.secondary[5]
  
  psi[6,1] <- (1-beta.stay.secondary[6])*(move[6,1]/(move[6,1]+move[6,2]+move[6,3]+move[6,4]+move[6,5]+move[6,7]+move[6,8]+move[6,9]+move[6,10]+move[6,11]+move[6,12]))
  psi[6,2] <- (1-beta.stay.secondary[6])*(move[6,2]/(move[6,1]+move[6,2]+move[6,3]+move[6,4]+move[6,5]+move[6,7]+move[6,8]+move[6,9]+move[6,10]+move[6,11]+move[6,12]))
  psi[6,3] <- (1-beta.stay.secondary[6])*(move[6,3]/(move[6,1]+move[6,2]+move[6,3]+move[6,4]+move[6,5]+move[6,7]+move[6,8]+move[6,9]+move[6,10]+move[6,11]+move[6,12]))
  psi[6,4] <- (1-beta.stay.secondary[6])*(move[6,4]/(move[6,1]+move[6,2]+move[6,3]+move[6,4]+move[6,5]+move[6,7]+move[6,8]+move[6,9]+move[6,10]+move[6,11]+move[6,12]))
  psi[6,5] <- (1-beta.stay.secondary[6])*(move[6,5]/(move[6,1]+move[6,2]+move[6,3]+move[6,4]+move[6,5]+move[6,7]+move[6,8]+move[6,9]+move[6,10]+move[6,11]+move[6,12]))
  psi[6,7] <- (1-beta.stay.secondary[6])*(move[6,7]/(move[6,1]+move[6,2]+move[6,3]+move[6,4]+move[6,5]+move[6,7]+move[6,8]+move[6,9]+move[6,10]+move[6,11]+move[6,12]))
  psi[6,8] <- (1-beta.stay.secondary[6])*(move[6,8]/(move[6,1]+move[6,2]+move[6,3]+move[6,4]+move[6,5]+move[6,7]+move[6,8]+move[6,9]+move[6,10]+move[6,11]+move[6,12]))
  psi[6,9] <- (1-beta.stay.secondary[6])*(move[6,9]/(move[6,1]+move[6,2]+move[6,3]+move[6,4]+move[6,5]+move[6,7]+move[6,8]+move[6,9]+move[6,10]+move[6,11]+move[6,12]))
  psi[6,10] <- (1-beta.stay.secondary[6])*(move[6,10]/(move[6,1]+move[6,2]+move[6,3]+move[6,4]+move[6,5]+move[6,7]+move[6,8]+move[6,9]+move[6,10]+move[6,11]+move[6,12]))
  psi[6,11] <- (1-beta.stay.secondary[6])*(move[6,11]/(move[6,1]+move[6,2]+move[6,3]+move[6,4]+move[6,5]+move[6,7]+move[6,8]+move[6,9]+move[6,10]+move[6,11]+move[6,12]))
  psi[6,12] <- (1-beta.stay.secondary[6])*(move[6,12]/(move[6,1]+move[6,2]+move[6,3]+move[6,4]+move[6,5]+move[6,7]+move[6,8]+move[6,9]+move[6,10]+move[6,11]+move[6,12]))
  psi[6,6] <- beta.stay.secondary[6]
  
  psi[7,1] <- (1-beta.stay.secondary[7])*(move[7,1]/(move[7,1]+move[7,2]+move[7,3]+move[7,4]+move[7,5]+move[7,6]+move[7,8]+move[7,9]+move[7,10]+move[7,11]+move[7,12]))
  psi[7,2] <- (1-beta.stay.secondary[7])*(move[7,2]/(move[7,1]+move[7,2]+move[7,3]+move[7,4]+move[7,5]+move[7,6]+move[7,8]+move[7,9]+move[7,10]+move[7,11]+move[7,12]))
  psi[7,3] <- (1-beta.stay.secondary[7])*(move[7,3]/(move[7,1]+move[7,2]+move[7,3]+move[7,4]+move[7,5]+move[7,6]+move[7,8]+move[7,9]+move[7,10]+move[7,11]+move[7,12]))
  psi[7,4] <- (1-beta.stay.secondary[7])*(move[7,4]/(move[7,1]+move[7,2]+move[7,3]+move[7,4]+move[7,5]+move[7,6]+move[7,8]+move[7,9]+move[7,10]+move[7,11]+move[7,12]))
  psi[7,5] <- (1-beta.stay.secondary[7])*(move[7,5]/(move[7,1]+move[7,2]+move[7,3]+move[7,4]+move[7,5]+move[7,6]+move[7,8]+move[7,9]+move[7,10]+move[7,11]+move[7,12]))
  psi[7,6] <- (1-beta.stay.secondary[7])*(move[7,6]/(move[7,1]+move[7,2]+move[7,3]+move[7,4]+move[7,5]+move[7,6]+move[7,8]+move[7,9]+move[7,10]+move[7,11]+move[7,12]))
  psi[7,8] <- (1-beta.stay.secondary[7])*(move[7,8]/(move[7,1]+move[7,2]+move[7,3]+move[7,4]+move[7,5]+move[7,6]+move[7,8]+move[7,9]+move[7,10]+move[7,11]+move[7,12]))
  psi[7,9] <- (1-beta.stay.secondary[7])*(move[7,9]/(move[7,1]+move[7,2]+move[7,3]+move[7,4]+move[7,5]+move[7,6]+move[7,8]+move[7,9]+move[7,10]+move[7,11]+move[7,12]))
  psi[7,10] <- (1-beta.stay.secondary[7])*(move[7,10]/(move[7,1]+move[7,2]+move[7,3]+move[7,4]+move[7,5]+move[7,6]+move[7,8]+move[7,9]+move[7,10]+move[7,11]+move[7,12]))
  psi[7,11] <- (1-beta.stay.secondary[7])*(move[7,11]/(move[7,1]+move[7,2]+move[7,3]+move[7,4]+move[7,5]+move[7,6]+move[7,8]+move[7,9]+move[7,10]+move[7,11]+move[7,12]))
  psi[7,12] <- (1-beta.stay.secondary[7])*(move[7,12]/(move[7,1]+move[7,2]+move[7,3]+move[7,4]+move[7,5]+move[7,6]+move[7,8]+move[7,9]+move[7,10]+move[7,11]+move[7,12]))
  psi[7,7] <- beta.stay.secondary[7]
  
  psi[8,1] <- (1-beta.stay.secondary[8])*(move[8,1]/(move[8,1]+move[8,2]+move[8,3]+move[8,4]+move[8,5]+move[8,6]+move[8,7]+move[8,9]+move[8,10]+move[8,11]+move[8,12]))
  psi[8,2] <- (1-beta.stay.secondary[8])*(move[8,2]/(move[8,1]+move[8,2]+move[8,3]+move[8,4]+move[8,5]+move[8,6]+move[8,7]+move[8,9]+move[8,10]+move[8,11]+move[8,12]))
  psi[8,3] <- (1-beta.stay.secondary[8])*(move[8,3]/(move[8,1]+move[8,2]+move[8,3]+move[8,4]+move[8,5]+move[8,6]+move[8,7]+move[8,9]+move[8,10]+move[8,11]+move[8,12]))
  psi[8,4] <- (1-beta.stay.secondary[8])*(move[8,4]/(move[8,1]+move[8,2]+move[8,3]+move[8,4]+move[8,5]+move[8,6]+move[8,7]+move[8,9]+move[8,10]+move[8,11]+move[8,12]))
  psi[8,5] <- (1-beta.stay.secondary[8])*(move[8,5]/(move[8,1]+move[8,2]+move[8,3]+move[8,4]+move[8,5]+move[8,6]+move[8,7]+move[8,9]+move[8,10]+move[8,11]+move[8,12]))
  psi[8,6] <- (1-beta.stay.secondary[8])*(move[8,6]/(move[8,1]+move[8,2]+move[8,3]+move[8,4]+move[8,5]+move[8,6]+move[8,7]+move[8,9]+move[8,10]+move[8,11]+move[8,12]))
  psi[8,7] <- (1-beta.stay.secondary[8])*(move[8,7]/(move[8,1]+move[8,2]+move[8,3]+move[8,4]+move[8,5]+move[8,6]+move[8,7]+move[8,9]+move[8,10]+move[8,11]+move[8,12]))
  psi[8,9] <- (1-beta.stay.secondary[8])*(move[8,9]/(move[8,1]+move[8,2]+move[8,3]+move[8,4]+move[8,5]+move[8,6]+move[8,7]+move[8,9]+move[8,10]+move[8,11]+move[8,12]))
  psi[8,10] <- (1-beta.stay.secondary[8])*(move[8,10]/(move[8,1]+move[8,2]+move[8,3]+move[8,4]+move[8,5]+move[8,6]+move[8,7]+move[8,9]+move[8,10]+move[8,11]+move[8,12]))
  psi[8,11] <- (1-beta.stay.secondary[8])*(move[8,11]/(move[8,1]+move[8,2]+move[8,3]+move[8,4]+move[8,5]+move[8,6]+move[8,7]+move[8,9]+move[8,10]+move[8,11]+move[8,12]))
  psi[8,12] <- (1-beta.stay.secondary[8])*(move[8,12]/(move[8,1]+move[8,2]+move[8,3]+move[8,4]+move[8,5]+move[8,6]+move[8,7]+move[8,9]+move[8,10]+move[8,11]+move[8,12]))
  psi[8,8] <- beta.stay.secondary[8]
  
  psi[9,1] <- (1-beta.stay.secondary[9])*(move[9,1]/(move[9,1]+move[9,2]+move[9,3]+move[9,4]+move[9,5]+move[9,6]+move[9,7]+move[9,8]+move[9,10]+move[9,11]+move[9,12]))
  psi[9,2] <- (1-beta.stay.secondary[9])*(move[9,2]/(move[9,1]+move[9,2]+move[9,3]+move[9,4]+move[9,5]+move[9,6]+move[9,7]+move[9,8]+move[9,10]+move[9,11]+move[9,12]))
  psi[9,3] <- (1-beta.stay.secondary[9])*(move[9,3]/(move[9,1]+move[9,2]+move[9,3]+move[9,4]+move[9,5]+move[9,6]+move[9,7]+move[9,8]+move[9,10]+move[9,11]+move[9,12]))
  psi[9,4] <- (1-beta.stay.secondary[9])*(move[9,4]/(move[9,1]+move[9,2]+move[9,3]+move[9,4]+move[9,5]+move[9,6]+move[9,7]+move[9,8]+move[9,10]+move[9,11]+move[9,12]))
  psi[9,5] <- (1-beta.stay.secondary[9])*(move[9,5]/(move[9,1]+move[9,2]+move[9,3]+move[9,4]+move[9,5]+move[9,6]+move[9,7]+move[9,8]+move[9,10]+move[9,11]+move[9,12]))
  psi[9,6] <- (1-beta.stay.secondary[9])*(move[9,6]/(move[9,1]+move[9,2]+move[9,3]+move[9,4]+move[9,5]+move[9,6]+move[9,7]+move[9,8]+move[9,10]+move[9,11]+move[9,12]))
  psi[9,7] <- (1-beta.stay.secondary[9])*(move[9,7]/(move[9,1]+move[9,2]+move[9,3]+move[9,4]+move[9,5]+move[9,6]+move[9,7]+move[9,8]+move[9,10]+move[9,11]+move[9,12]))
  psi[9,8] <- (1-beta.stay.secondary[9])*(move[9,8]/(move[9,1]+move[9,2]+move[9,3]+move[9,4]+move[9,5]+move[9,6]+move[9,7]+move[9,8]+move[9,10]+move[9,11]+move[9,12]))
  psi[9,10] <- (1-beta.stay.secondary[9])*(move[9,10]/(move[9,1]+move[9,2]+move[9,3]+move[9,4]+move[9,5]+move[9,6]+move[9,7]+move[9,8]+move[9,10]+move[9,11]+move[9,12]))
  psi[9,11] <- (1-beta.stay.secondary[9])*(move[9,11]/(move[9,1]+move[9,2]+move[9,3]+move[9,4]+move[9,5]+move[9,6]+move[9,7]+move[9,8]+move[9,10]+move[9,11]+move[9,12]))
  psi[9,12] <- (1-beta.stay.secondary[9])*(move[9,12]/(move[9,1]+move[9,2]+move[9,3]+move[9,4]+move[9,5]+move[9,6]+move[9,7]+move[9,8]+move[9,10]+move[9,11]+move[9,12]))
  psi[9,9] <- beta.stay.secondary[9]
  
  psi[10,1] <- (1-beta.stay.secondary[10])*(move[10,1]/(move[10,1]+move[10,2]+move[10,3]+move[10,4]+move[10,5]+move[10,6]+move[10,7]+move[10,8]+move[10,9]+move[10,11]+move[10,12]))
  psi[10,2] <- (1-beta.stay.secondary[10])*(move[10,2]/(move[10,1]+move[10,2]+move[10,3]+move[10,4]+move[10,5]+move[10,6]+move[10,7]+move[10,8]+move[10,9]+move[10,11]+move[10,12]))
  psi[10,3] <- (1-beta.stay.secondary[10])*(move[10,3]/(move[10,1]+move[10,2]+move[10,3]+move[10,4]+move[10,5]+move[10,6]+move[10,7]+move[10,8]+move[10,9]+move[10,11]+move[10,12]))
  psi[10,4] <- (1-beta.stay.secondary[10])*(move[10,4]/(move[10,1]+move[10,2]+move[10,3]+move[10,4]+move[10,5]+move[10,6]+move[10,7]+move[10,8]+move[10,9]+move[10,11]+move[10,12]))
  psi[10,5] <- (1-beta.stay.secondary[10])*(move[10,5]/(move[10,1]+move[10,2]+move[10,3]+move[10,4]+move[10,5]+move[10,6]+move[10,7]+move[10,8]+move[10,9]+move[10,11]+move[10,12]))
  psi[10,6] <- (1-beta.stay.secondary[10])*(move[10,6]/(move[10,1]+move[10,2]+move[10,3]+move[10,4]+move[10,5]+move[10,6]+move[10,7]+move[10,8]+move[10,9]+move[10,11]+move[10,12]))
  psi[10,7] <- (1-beta.stay.secondary[10])*(move[10,7]/(move[10,1]+move[10,2]+move[10,3]+move[10,4]+move[10,5]+move[10,6]+move[10,7]+move[10,8]+move[10,9]+move[10,11]+move[10,12]))
  psi[10,8] <- (1-beta.stay.secondary[10])*(move[10,8]/(move[10,1]+move[10,2]+move[10,3]+move[10,4]+move[10,5]+move[10,6]+move[10,7]+move[10,8]+move[10,9]+move[10,11]+move[10,12]))
  psi[10,9] <- (1-beta.stay.secondary[10])*(move[10,9]/(move[10,1]+move[10,2]+move[10,3]+move[10,4]+move[10,5]+move[10,6]+move[10,7]+move[10,8]+move[10,9]+move[10,11]+move[10,12]))
  psi[10,11] <- (1-beta.stay.secondary[10])*(move[10,11]/(move[10,1]+move[10,2]+move[10,3]+move[10,4]+move[10,5]+move[10,6]+move[10,7]+move[10,8]+move[10,9]+move[10,11]+move[10,12]))
  psi[10,12] <- (1-beta.stay.secondary[10])*(move[10,12]/(move[10,1]+move[10,2]+move[10,3]+move[10,4]+move[10,5]+move[10,6]+move[10,7]+move[10,8]+move[10,9]+move[10,11]+move[10,12]))
  psi[10,10] <- beta.stay.secondary[10]
  
  psi[11,1] <- (1-beta.stay.secondary[11])*(move[11,1]/(move[11,1]+move[11,2]+move[11,3]+move[11,4]+move[11,5]+move[11,6]+move[11,7]+move[11,8]+move[11,9]+move[11,10]+move[11,12]))
  psi[11,2] <- (1-beta.stay.secondary[11])*(move[11,2]/(move[11,1]+move[11,2]+move[11,3]+move[11,4]+move[11,5]+move[11,6]+move[11,7]+move[11,8]+move[11,9]+move[11,10]+move[11,12]))
  psi[11,3] <- (1-beta.stay.secondary[11])*(move[11,3]/(move[11,1]+move[11,2]+move[11,3]+move[11,4]+move[11,5]+move[11,6]+move[11,7]+move[11,8]+move[11,9]+move[11,10]+move[11,12]))
  psi[11,4] <- (1-beta.stay.secondary[11])*(move[11,4]/(move[11,1]+move[11,2]+move[11,3]+move[11,4]+move[11,5]+move[11,6]+move[11,7]+move[11,8]+move[11,9]+move[11,10]+move[11,12]))
  psi[11,5] <- (1-beta.stay.secondary[11])*(move[11,5]/(move[11,1]+move[11,2]+move[11,3]+move[11,4]+move[11,5]+move[11,6]+move[11,7]+move[11,8]+move[11,9]+move[11,10]+move[11,12]))
  psi[11,6] <- (1-beta.stay.secondary[11])*(move[11,6]/(move[11,1]+move[11,2]+move[11,3]+move[11,4]+move[11,5]+move[11,6]+move[11,7]+move[11,8]+move[11,9]+move[11,10]+move[11,12]))
  psi[11,7] <- (1-beta.stay.secondary[11])*(move[11,7]/(move[11,1]+move[11,2]+move[11,3]+move[11,4]+move[11,5]+move[11,6]+move[11,7]+move[11,8]+move[11,9]+move[11,10]+move[11,12]))
  psi[11,8] <- (1-beta.stay.secondary[11])*(move[11,8]/(move[11,1]+move[11,2]+move[11,3]+move[11,4]+move[11,5]+move[11,6]+move[11,7]+move[11,8]+move[11,9]+move[11,10]+move[11,12]))
  psi[11,9] <- (1-beta.stay.secondary[11])*(move[11,9]/(move[11,1]+move[11,2]+move[11,3]+move[11,4]+move[11,5]+move[11,6]+move[11,7]+move[11,8]+move[11,9]+move[11,10]+move[11,12]))
  psi[11,10] <- (1-beta.stay.secondary[11])*(move[11,10]/(move[11,1]+move[11,2]+move[11,3]+move[11,4]+move[11,5]+move[11,6]+move[11,7]+move[11,8]+move[11,9]+move[11,10]+move[11,12]))
  psi[11,12] <- (1-beta.stay.secondary[11])*(move[11,12]/(move[11,1]+move[11,2]+move[11,3]+move[11,4]+move[11,5]+move[11,6]+move[11,7]+move[11,8]+move[11,9]+move[11,10]+move[11,12]))
  psi[11,11] <- beta.stay.secondary[11]
  
  psi[12,1] <- (1-beta.stay.secondary[12])*(move[12,1]/(move[12,1]+move[12,2]+move[12,3]+move[12,4]+move[12,5]+move[12,6]+move[12,7]+move[12,8]+move[12,9]+move[12,10]+move[12,11]))
  psi[12,2] <- (1-beta.stay.secondary[12])*(move[12,2]/(move[12,1]+move[12,2]+move[12,3]+move[12,4]+move[12,5]+move[12,6]+move[12,7]+move[12,8]+move[12,9]+move[12,10]+move[12,11]))
  psi[12,3] <- (1-beta.stay.secondary[12])*(move[12,3]/(move[12,1]+move[12,2]+move[12,3]+move[12,4]+move[12,5]+move[12,6]+move[12,7]+move[12,8]+move[12,9]+move[12,10]+move[12,11]))
  psi[12,4] <- (1-beta.stay.secondary[12])*(move[12,4]/(move[12,1]+move[12,2]+move[12,3]+move[12,4]+move[12,5]+move[12,6]+move[12,7]+move[12,8]+move[12,9]+move[12,10]+move[12,11]))
  psi[12,5] <- (1-beta.stay.secondary[12])*(move[12,5]/(move[12,1]+move[12,2]+move[12,3]+move[12,4]+move[12,5]+move[12,6]+move[12,7]+move[12,8]+move[12,9]+move[12,10]+move[12,11]))
  psi[12,6] <- (1-beta.stay.secondary[12])*(move[12,6]/(move[12,1]+move[12,2]+move[12,3]+move[12,4]+move[12,5]+move[12,6]+move[12,7]+move[12,8]+move[12,9]+move[12,10]+move[12,11]))
  psi[12,7] <- (1-beta.stay.secondary[12])*(move[12,7]/(move[12,1]+move[12,2]+move[12,3]+move[12,4]+move[12,5]+move[12,6]+move[12,7]+move[12,8]+move[12,9]+move[12,10]+move[12,11]))
  psi[12,8] <- (1-beta.stay.secondary[12])*(move[12,8]/(move[12,1]+move[12,2]+move[12,3]+move[12,4]+move[12,5]+move[12,6]+move[12,7]+move[12,8]+move[12,9]+move[12,10]+move[12,11]))
  psi[12,9] <- (1-beta.stay.secondary[12])*(move[12,9]/(move[12,1]+move[12,2]+move[12,3]+move[12,4]+move[12,5]+move[12,6]+move[12,7]+move[12,8]+move[12,9]+move[12,10]+move[12,11]))
  psi[12,10] <- (1-beta.stay.secondary[12])*(move[12,10]/(move[12,1]+move[12,2]+move[12,3]+move[12,4]+move[12,5]+move[12,6]+move[12,7]+move[12,8]+move[12,9]+move[12,10]+move[12,11]))
  psi[12,11] <- (1-beta.stay.secondary[12])*(move[12,11]/(move[12,1]+move[12,2]+move[12,3]+move[12,4]+move[12,5]+move[12,6]+move[12,7]+move[12,8]+move[12,9]+move[12,10]+move[12,11]))
  psi[12,12] <- beta.stay.secondary[12]
  
  # retention probabilities
  phiintercept <- rep(param[45],T)
  phigradient <- -exp(param[46:51])
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
  sigma.primary <- exp(param[52])

  Fid <- matrix(0, 12, 1)
  Fid <- 1/(1+exp(-param[53:64]))
  
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
  fid[1,2] <- (1-Fid[1])*(move.primary[1,2]/(move.primary[1,2]+move.primary[1,3]+move.primary[1,4]+move.primary[1,5]+move.primary[1,6]+move.primary[1,7]+move.primary[1,8]+move.primary[1,9]+move.primary[1,10]+move.primary[1,11]+move.primary[1,12]))
  fid[1,3] <- (1-Fid[1])*(move.primary[1,3]/(move.primary[1,2]+move.primary[1,3]+move.primary[1,4]+move.primary[1,5]+move.primary[1,6]+move.primary[1,7]+move.primary[1,8]+move.primary[1,9]+move.primary[1,10]+move.primary[1,11]+move.primary[1,12]))
  fid[1,4] <- (1-Fid[1])*(move.primary[1,4]/(move.primary[1,2]+move.primary[1,3]+move.primary[1,4]+move.primary[1,5]+move.primary[1,6]+move.primary[1,7]+move.primary[1,8]+move.primary[1,9]+move.primary[1,10]+move.primary[1,11]+move.primary[1,12]))
  fid[1,5] <- (1-Fid[1])*(move.primary[1,5]/(move.primary[1,2]+move.primary[1,3]+move.primary[1,4]+move.primary[1,5]+move.primary[1,6]+move.primary[1,7]+move.primary[1,8]+move.primary[1,9]+move.primary[1,10]+move.primary[1,11]+move.primary[1,12]))
  fid[1,6] <- (1-Fid[1])*(move.primary[1,6]/(move.primary[1,2]+move.primary[1,3]+move.primary[1,4]+move.primary[1,5]+move.primary[1,6]+move.primary[1,7]+move.primary[1,8]+move.primary[1,9]+move.primary[1,10]+move.primary[1,11]+move.primary[1,12]))
  fid[1,7] <- (1-Fid[1])*(move.primary[1,7]/(move.primary[1,2]+move.primary[1,3]+move.primary[1,4]+move.primary[1,5]+move.primary[1,6]+move.primary[1,7]+move.primary[1,8]+move.primary[1,9]+move.primary[1,10]+move.primary[1,11]+move.primary[1,12]))
  fid[1,8] <- (1-Fid[1])*(move.primary[1,8]/(move.primary[1,2]+move.primary[1,3]+move.primary[1,4]+move.primary[1,5]+move.primary[1,6]+move.primary[1,7]+move.primary[1,8]+move.primary[1,9]+move.primary[1,10]+move.primary[1,11]+move.primary[1,12]))
  fid[1,9] <- (1-Fid[1])*(move.primary[1,9]/(move.primary[1,2]+move.primary[1,3]+move.primary[1,4]+move.primary[1,5]+move.primary[1,6]+move.primary[1,7]+move.primary[1,8]+move.primary[1,9]+move.primary[1,10]+move.primary[1,11]+move.primary[1,12]))
  fid[1,10] <- (1-Fid[1])*(move.primary[1,10]/(move.primary[1,2]+move.primary[1,3]+move.primary[1,4]+move.primary[1,5]+move.primary[1,6]+move.primary[1,7]+move.primary[1,8]+move.primary[1,9]+move.primary[1,10]+move.primary[1,11]+move.primary[1,12]))
  fid[1,11] <- (1-Fid[1])*(move.primary[1,11]/(move.primary[1,2]+move.primary[1,3]+move.primary[1,4]+move.primary[1,5]+move.primary[1,6]+move.primary[1,7]+move.primary[1,8]+move.primary[1,9]+move.primary[1,10]+move.primary[1,11]+move.primary[1,12]))
  fid[1,12] <- (1-Fid[1])*(move.primary[1,12]/(move.primary[1,2]+move.primary[1,3]+move.primary[1,4]+move.primary[1,5]+move.primary[1,6]+move.primary[1,7]+move.primary[1,8]+move.primary[1,9]+move.primary[1,10]+move.primary[1,11]+move.primary[1,12]))
  fid[1,1] <- Fid[1]
  
  fid[2,1] <- (1-Fid[2])*(move.primary[2,1]/(move.primary[2,1]+move.primary[2,3]+move.primary[2,4]+move.primary[2,5]+move.primary[2,6]+move.primary[2,7]+move.primary[2,8]+move.primary[2,9]+move.primary[2,10]+move.primary[2,11]+move.primary[2,12]))
  fid[2,3] <- (1-Fid[2])*(move.primary[2,3]/(move.primary[2,1]+move.primary[2,3]+move.primary[2,4]+move.primary[2,5]+move.primary[2,6]+move.primary[2,7]+move.primary[2,8]+move.primary[2,9]+move.primary[2,10]+move.primary[2,11]+move.primary[2,12]))
  fid[2,4] <- (1-Fid[2])*(move.primary[2,4]/(move.primary[2,1]+move.primary[2,3]+move.primary[2,4]+move.primary[2,5]+move.primary[2,6]+move.primary[2,7]+move.primary[2,8]+move.primary[2,9]+move.primary[2,10]+move.primary[2,11]+move.primary[2,12]))
  fid[2,5] <- (1-Fid[2])*(move.primary[2,5]/(move.primary[2,1]+move.primary[2,3]+move.primary[2,4]+move.primary[2,5]+move.primary[2,6]+move.primary[2,7]+move.primary[2,8]+move.primary[2,9]+move.primary[2,10]+move.primary[2,11]+move.primary[2,12]))
  fid[2,6] <- (1-Fid[2])*(move.primary[2,6]/(move.primary[2,1]+move.primary[2,3]+move.primary[2,4]+move.primary[2,5]+move.primary[2,6]+move.primary[2,7]+move.primary[2,8]+move.primary[2,9]+move.primary[2,10]+move.primary[2,11]+move.primary[2,12]))
  fid[2,7] <- (1-Fid[2])*(move.primary[2,7]/(move.primary[2,1]+move.primary[2,3]+move.primary[2,4]+move.primary[2,5]+move.primary[2,6]+move.primary[2,7]+move.primary[2,8]+move.primary[2,9]+move.primary[2,10]+move.primary[2,11]+move.primary[2,12]))
  fid[2,8] <- (1-Fid[2])*(move.primary[2,8]/(move.primary[2,1]+move.primary[2,3]+move.primary[2,4]+move.primary[2,5]+move.primary[2,6]+move.primary[2,7]+move.primary[2,8]+move.primary[2,9]+move.primary[2,10]+move.primary[2,11]+move.primary[2,12]))
  fid[2,9] <- (1-Fid[2])*(move.primary[2,9]/(move.primary[2,1]+move.primary[2,3]+move.primary[2,4]+move.primary[2,5]+move.primary[2,6]+move.primary[2,7]+move.primary[2,8]+move.primary[2,9]+move.primary[2,10]+move.primary[2,11]+move.primary[2,12]))
  fid[2,10] <- (1-Fid[2])*(move.primary[2,10]/(move.primary[2,1]+move.primary[2,3]+move.primary[2,4]+move.primary[2,5]+move.primary[2,6]+move.primary[2,7]+move.primary[2,8]+move.primary[2,9]+move.primary[2,10]+move.primary[2,11]+move.primary[2,12]))
  fid[2,11] <- (1-Fid[2])*(move.primary[2,11]/(move.primary[2,1]+move.primary[2,3]+move.primary[2,4]+move.primary[2,5]+move.primary[2,6]+move.primary[2,7]+move.primary[2,8]+move.primary[2,9]+move.primary[2,10]+move.primary[2,11]+move.primary[2,12]))
  fid[2,12] <- (1-Fid[2])*(move.primary[2,12]/(move.primary[2,1]+move.primary[2,3]+move.primary[2,4]+move.primary[2,5]+move.primary[2,6]+move.primary[2,7]+move.primary[2,8]+move.primary[2,9]+move.primary[2,10]+move.primary[2,11]+move.primary[2,12]))
  fid[2,2] <- Fid[2]
  
  fid[3,1] <- (1-Fid[3])*(move.primary[3,1]/(move.primary[3,1]+move.primary[3,2]+move.primary[3,4]+move.primary[3,5]+move.primary[3,6]+move.primary[3,7]+move.primary[3,8]+move.primary[3,9]+move.primary[3,10]+move.primary[3,11]+move.primary[3,12]))
  fid[3,2] <- (1-Fid[3])*(move.primary[3,2]/(move.primary[3,1]+move.primary[3,2]+move.primary[3,4]+move.primary[3,5]+move.primary[3,6]+move.primary[3,7]+move.primary[3,8]+move.primary[3,9]+move.primary[3,10]+move.primary[3,11]+move.primary[3,12]))
  fid[3,4] <- (1-Fid[3])*(move.primary[3,4]/(move.primary[3,1]+move.primary[3,2]+move.primary[3,4]+move.primary[3,5]+move.primary[3,6]+move.primary[3,7]+move.primary[3,8]+move.primary[3,9]+move.primary[3,10]+move.primary[3,11]+move.primary[3,12]))
  fid[3,5] <- (1-Fid[3])*(move.primary[3,5]/(move.primary[3,1]+move.primary[3,2]+move.primary[3,4]+move.primary[3,5]+move.primary[3,6]+move.primary[3,7]+move.primary[3,8]+move.primary[3,9]+move.primary[3,10]+move.primary[3,11]+move.primary[3,12]))
  fid[3,6] <- (1-Fid[3])*(move.primary[3,6]/(move.primary[3,1]+move.primary[3,2]+move.primary[3,4]+move.primary[3,5]+move.primary[3,6]+move.primary[3,7]+move.primary[3,8]+move.primary[3,9]+move.primary[3,10]+move.primary[3,11]+move.primary[3,12]))
  fid[3,7] <- (1-Fid[3])*(move.primary[3,7]/(move.primary[3,1]+move.primary[3,2]+move.primary[3,4]+move.primary[3,5]+move.primary[3,6]+move.primary[3,7]+move.primary[3,8]+move.primary[3,9]+move.primary[3,10]+move.primary[3,11]+move.primary[3,12]))
  fid[3,8] <- (1-Fid[3])*(move.primary[3,8]/(move.primary[3,1]+move.primary[3,2]+move.primary[3,4]+move.primary[3,5]+move.primary[3,6]+move.primary[3,7]+move.primary[3,8]+move.primary[3,9]+move.primary[3,10]+move.primary[3,11]+move.primary[3,12]))
  fid[3,9] <- (1-Fid[3])*(move.primary[3,9]/(move.primary[3,1]+move.primary[3,2]+move.primary[3,4]+move.primary[3,5]+move.primary[3,6]+move.primary[3,7]+move.primary[3,8]+move.primary[3,9]+move.primary[3,10]+move.primary[3,11]+move.primary[3,12]))
  fid[3,10] <- (1-Fid[3])*(move.primary[3,10]/(move.primary[3,1]+move.primary[3,2]+move.primary[3,4]+move.primary[3,5]+move.primary[3,6]+move.primary[3,7]+move.primary[3,8]+move.primary[3,9]+move.primary[3,10]+move.primary[3,11]+move.primary[3,12]))
  fid[3,11] <- (1-Fid[3])*(move.primary[3,11]/(move.primary[3,1]+move.primary[3,2]+move.primary[3,4]+move.primary[3,5]+move.primary[3,6]+move.primary[3,7]+move.primary[3,8]+move.primary[3,9]+move.primary[3,10]+move.primary[3,11]+move.primary[3,12]))
  fid[3,12] <- (1-Fid[3])*(move.primary[3,12]/(move.primary[3,1]+move.primary[3,2]+move.primary[3,4]+move.primary[3,5]+move.primary[3,6]+move.primary[3,7]+move.primary[3,8]+move.primary[3,9]+move.primary[3,10]+move.primary[3,11]+move.primary[3,12]))
  fid[3,3] <- Fid[3]
  
  fid[4,1] <- (1-Fid[4])*(move.primary[4,1]/(move.primary[4,1]+move.primary[4,2]+move.primary[4,3]+move.primary[4,5]+move.primary[4,6]+move.primary[4,7]+move.primary[4,8]+move.primary[4,9]+move.primary[4,10]+move.primary[4,11]+move.primary[4,12]))
  fid[4,2] <- (1-Fid[4])*(move.primary[4,2]/(move.primary[4,1]+move.primary[4,2]+move.primary[4,3]+move.primary[4,5]+move.primary[4,6]+move.primary[4,7]+move.primary[4,8]+move.primary[4,9]+move.primary[4,10]+move.primary[4,11]+move.primary[4,12]))
  fid[4,3] <- (1-Fid[4])*(move.primary[4,3]/(move.primary[4,1]+move.primary[4,2]+move.primary[4,3]+move.primary[4,5]+move.primary[4,6]+move.primary[4,7]+move.primary[4,8]+move.primary[4,9]+move.primary[4,10]+move.primary[4,11]+move.primary[4,12]))
  fid[4,5] <- (1-Fid[4])*(move.primary[4,5]/(move.primary[4,1]+move.primary[4,2]+move.primary[4,3]+move.primary[4,5]+move.primary[4,6]+move.primary[4,7]+move.primary[4,8]+move.primary[4,9]+move.primary[4,10]+move.primary[4,11]+move.primary[4,12]))
  fid[4,6] <- (1-Fid[4])*(move.primary[4,6]/(move.primary[4,1]+move.primary[4,2]+move.primary[4,3]+move.primary[4,5]+move.primary[4,6]+move.primary[4,7]+move.primary[4,8]+move.primary[4,9]+move.primary[4,10]+move.primary[4,11]+move.primary[4,12]))
  fid[4,7] <- (1-Fid[4])*(move.primary[4,7]/(move.primary[4,1]+move.primary[4,2]+move.primary[4,3]+move.primary[4,5]+move.primary[4,6]+move.primary[4,7]+move.primary[4,8]+move.primary[4,9]+move.primary[4,10]+move.primary[4,11]+move.primary[4,12]))
  fid[4,8] <- (1-Fid[4])*(move.primary[4,8]/(move.primary[4,1]+move.primary[4,2]+move.primary[4,3]+move.primary[4,5]+move.primary[4,6]+move.primary[4,7]+move.primary[4,8]+move.primary[4,9]+move.primary[4,10]+move.primary[4,11]+move.primary[4,12]))
  fid[4,9] <- (1-Fid[4])*(move.primary[4,9]/(move.primary[4,1]+move.primary[4,2]+move.primary[4,3]+move.primary[4,5]+move.primary[4,6]+move.primary[4,7]+move.primary[4,8]+move.primary[4,9]+move.primary[4,10]+move.primary[4,11]+move.primary[4,12]))
  fid[4,10] <- (1-Fid[4])*(move.primary[4,10]/(move.primary[4,1]+move.primary[4,2]+move.primary[4,3]+move.primary[4,5]+move.primary[4,6]+move.primary[4,7]+move.primary[4,8]+move.primary[4,9]+move.primary[4,10]+move.primary[4,11]+move.primary[4,12]))
  fid[4,11] <- (1-Fid[4])*(move.primary[4,11]/(move.primary[4,1]+move.primary[4,2]+move.primary[4,3]+move.primary[4,5]+move.primary[4,6]+move.primary[4,7]+move.primary[4,8]+move.primary[4,9]+move.primary[4,10]+move.primary[4,11]+move.primary[4,12]))
  fid[4,12] <- (1-Fid[4])*(move.primary[4,12]/(move.primary[4,1]+move.primary[4,2]+move.primary[4,3]+move.primary[4,5]+move.primary[4,6]+move.primary[4,7]+move.primary[4,8]+move.primary[4,9]+move.primary[4,10]+move.primary[4,11]+move.primary[4,12]))
  fid[4,4] <- Fid[4]
  
  fid[5,1] <- (1-Fid[5])*(move.primary[5,1]/(move.primary[5,1]+move.primary[5,2]+move.primary[5,3]+move.primary[5,4]+move.primary[5,6]+move.primary[5,7]+move.primary[5,8]+move.primary[5,9]+move.primary[5,10]+move.primary[5,11]+move.primary[5,12]))
  fid[5,2] <- (1-Fid[5])*(move.primary[5,2]/(move.primary[5,1]+move.primary[5,2]+move.primary[5,3]+move.primary[5,4]+move.primary[5,6]+move.primary[5,7]+move.primary[5,8]+move.primary[5,9]+move.primary[5,10]+move.primary[5,11]+move.primary[5,12]))
  fid[5,3] <- (1-Fid[5])*(move.primary[5,3]/(move.primary[5,1]+move.primary[5,2]+move.primary[5,3]+move.primary[5,4]+move.primary[5,6]+move.primary[5,7]+move.primary[5,8]+move.primary[5,9]+move.primary[5,10]+move.primary[5,11]+move.primary[5,12]))
  fid[5,4] <- (1-Fid[5])*(move.primary[5,4]/(move.primary[5,1]+move.primary[5,2]+move.primary[5,3]+move.primary[5,4]+move.primary[5,6]+move.primary[5,7]+move.primary[5,8]+move.primary[5,9]+move.primary[5,10]+move.primary[5,11]+move.primary[5,12]))
  fid[5,6] <- (1-Fid[5])*(move.primary[5,6]/(move.primary[5,1]+move.primary[5,2]+move.primary[5,3]+move.primary[5,4]+move.primary[5,6]+move.primary[5,7]+move.primary[5,8]+move.primary[5,9]+move.primary[5,10]+move.primary[5,11]+move.primary[5,12]))
  fid[5,7] <- (1-Fid[5])*(move.primary[5,7]/(move.primary[5,1]+move.primary[5,2]+move.primary[5,3]+move.primary[5,4]+move.primary[5,6]+move.primary[5,7]+move.primary[5,8]+move.primary[5,9]+move.primary[5,10]+move.primary[5,11]+move.primary[5,12]))
  fid[5,8] <- (1-Fid[5])*(move.primary[5,8]/(move.primary[5,1]+move.primary[5,2]+move.primary[5,3]+move.primary[5,4]+move.primary[5,6]+move.primary[5,7]+move.primary[5,8]+move.primary[5,9]+move.primary[5,10]+move.primary[5,11]+move.primary[5,12]))
  fid[5,9] <- (1-Fid[5])*(move.primary[5,9]/(move.primary[5,1]+move.primary[5,2]+move.primary[5,3]+move.primary[5,4]+move.primary[5,6]+move.primary[5,7]+move.primary[5,8]+move.primary[5,9]+move.primary[5,10]+move.primary[5,11]+move.primary[5,12]))
  fid[5,10] <- (1-Fid[5])*(move.primary[5,10]/(move.primary[5,1]+move.primary[5,2]+move.primary[5,3]+move.primary[5,4]+move.primary[5,6]+move.primary[5,7]+move.primary[5,8]+move.primary[5,9]+move.primary[5,10]+move.primary[5,11]+move.primary[5,12]))
  fid[5,11] <- (1-Fid[5])*(move.primary[5,11]/(move.primary[5,1]+move.primary[5,2]+move.primary[5,3]+move.primary[5,4]+move.primary[5,6]+move.primary[5,7]+move.primary[5,8]+move.primary[5,9]+move.primary[5,10]+move.primary[5,11]+move.primary[5,12]))
  fid[5,12] <- (1-Fid[5])*(move.primary[5,12]/(move.primary[5,1]+move.primary[5,2]+move.primary[5,3]+move.primary[5,4]+move.primary[5,6]+move.primary[5,7]+move.primary[5,8]+move.primary[5,9]+move.primary[5,10]+move.primary[5,11]+move.primary[5,12]))
  fid[5,5] <- Fid[5]
  
  fid[6,1] <- (1-Fid[6])*(move.primary[6,1]/(move.primary[6,1]+move.primary[6,2]+move.primary[6,3]+move.primary[6,4]+move.primary[6,5]+move.primary[6,7]+move.primary[6,8]+move.primary[6,9]+move.primary[6,10]+move.primary[6,11]+move.primary[6,12]))
  fid[6,2] <- (1-Fid[6])*(move.primary[6,2]/(move.primary[6,1]+move.primary[6,2]+move.primary[6,3]+move.primary[6,4]+move.primary[6,5]+move.primary[6,7]+move.primary[6,8]+move.primary[6,9]+move.primary[6,10]+move.primary[6,11]+move.primary[6,12]))
  fid[6,3] <- (1-Fid[6])*(move.primary[6,3]/(move.primary[6,1]+move.primary[6,2]+move.primary[6,3]+move.primary[6,4]+move.primary[6,5]+move.primary[6,7]+move.primary[6,8]+move.primary[6,9]+move.primary[6,10]+move.primary[6,11]+move.primary[6,12]))
  fid[6,4] <- (1-Fid[6])*(move.primary[6,4]/(move.primary[6,1]+move.primary[6,2]+move.primary[6,3]+move.primary[6,4]+move.primary[6,5]+move.primary[6,7]+move.primary[6,8]+move.primary[6,9]+move.primary[6,10]+move.primary[6,11]+move.primary[6,12]))
  fid[6,5] <- (1-Fid[6])*(move.primary[6,5]/(move.primary[6,1]+move.primary[6,2]+move.primary[6,3]+move.primary[6,4]+move.primary[6,5]+move.primary[6,7]+move.primary[6,8]+move.primary[6,9]+move.primary[6,10]+move.primary[6,11]+move.primary[6,12]))
  fid[6,7] <- (1-Fid[6])*(move.primary[6,7]/(move.primary[6,1]+move.primary[6,2]+move.primary[6,3]+move.primary[6,4]+move.primary[6,5]+move.primary[6,7]+move.primary[6,8]+move.primary[6,9]+move.primary[6,10]+move.primary[6,11]+move.primary[6,12]))
  fid[6,8] <- (1-Fid[6])*(move.primary[6,8]/(move.primary[6,1]+move.primary[6,2]+move.primary[6,3]+move.primary[6,4]+move.primary[6,5]+move.primary[6,7]+move.primary[6,8]+move.primary[6,9]+move.primary[6,10]+move.primary[6,11]+move.primary[6,12]))
  fid[6,9] <- (1-Fid[6])*(move.primary[6,9]/(move.primary[6,1]+move.primary[6,2]+move.primary[6,3]+move.primary[6,4]+move.primary[6,5]+move.primary[6,7]+move.primary[6,8]+move.primary[6,9]+move.primary[6,10]+move.primary[6,11]+move.primary[6,12]))
  fid[6,10] <- (1-Fid[6])*(move.primary[6,10]/(move.primary[6,1]+move.primary[6,2]+move.primary[6,3]+move.primary[6,4]+move.primary[6,5]+move.primary[6,7]+move.primary[6,8]+move.primary[6,9]+move.primary[6,10]+move.primary[6,11]+move.primary[6,12]))
  fid[6,11] <- (1-Fid[6])*(move.primary[6,11]/(move.primary[6,1]+move.primary[6,2]+move.primary[6,3]+move.primary[6,4]+move.primary[6,5]+move.primary[6,7]+move.primary[6,8]+move.primary[6,9]+move.primary[6,10]+move.primary[6,11]+move.primary[6,12]))
  fid[6,12] <- (1-Fid[6])*(move.primary[6,12]/(move.primary[6,1]+move.primary[6,2]+move.primary[6,3]+move.primary[6,4]+move.primary[6,5]+move.primary[6,7]+move.primary[6,8]+move.primary[6,9]+move.primary[6,10]+move.primary[6,11]+move.primary[6,12]))
  fid[6,6] <- Fid[6]
  
  fid[7,1] <- (1-Fid[7])*(move.primary[7,1]/(move.primary[7,1]+move.primary[7,2]+move.primary[7,3]+move.primary[7,4]+move.primary[7,5]+move.primary[7,6]+move.primary[7,8]+move.primary[7,9]+move.primary[7,10]+move.primary[7,11]+move.primary[7,12]))
  fid[7,2] <- (1-Fid[7])*(move.primary[7,2]/(move.primary[7,1]+move.primary[7,2]+move.primary[7,3]+move.primary[7,4]+move.primary[7,5]+move.primary[7,6]+move.primary[7,8]+move.primary[7,9]+move.primary[7,10]+move.primary[7,11]+move.primary[7,12]))
  fid[7,3] <- (1-Fid[7])*(move.primary[7,3]/(move.primary[7,1]+move.primary[7,2]+move.primary[7,3]+move.primary[7,4]+move.primary[7,5]+move.primary[7,6]+move.primary[7,8]+move.primary[7,9]+move.primary[7,10]+move.primary[7,11]+move.primary[7,12]))
  fid[7,4] <- (1-Fid[7])*(move.primary[7,4]/(move.primary[7,1]+move.primary[7,2]+move.primary[7,3]+move.primary[7,4]+move.primary[7,5]+move.primary[7,6]+move.primary[7,8]+move.primary[7,9]+move.primary[7,10]+move.primary[7,11]+move.primary[7,12]))
  fid[7,5] <- (1-Fid[7])*(move.primary[7,5]/(move.primary[7,1]+move.primary[7,2]+move.primary[7,3]+move.primary[7,4]+move.primary[7,5]+move.primary[7,6]+move.primary[7,8]+move.primary[7,9]+move.primary[7,10]+move.primary[7,11]+move.primary[7,12]))
  fid[7,6] <- (1-Fid[7])*(move.primary[7,6]/(move.primary[7,1]+move.primary[7,2]+move.primary[7,3]+move.primary[7,4]+move.primary[7,5]+move.primary[7,6]+move.primary[7,8]+move.primary[7,9]+move.primary[7,10]+move.primary[7,11]+move.primary[7,12]))
  fid[7,8] <- (1-Fid[7])*(move.primary[7,8]/(move.primary[7,1]+move.primary[7,2]+move.primary[7,3]+move.primary[7,4]+move.primary[7,5]+move.primary[7,6]+move.primary[7,8]+move.primary[7,9]+move.primary[7,10]+move.primary[7,11]+move.primary[7,12]))
  fid[7,9] <- (1-Fid[7])*(move.primary[7,9]/(move.primary[7,1]+move.primary[7,2]+move.primary[7,3]+move.primary[7,4]+move.primary[7,5]+move.primary[7,6]+move.primary[7,8]+move.primary[7,9]+move.primary[7,10]+move.primary[7,11]+move.primary[7,12]))
  fid[7,10] <- (1-Fid[7])*(move.primary[7,10]/(move.primary[7,1]+move.primary[7,2]+move.primary[7,3]+move.primary[7,4]+move.primary[7,5]+move.primary[7,6]+move.primary[7,8]+move.primary[7,9]+move.primary[7,10]+move.primary[7,11]+move.primary[7,12]))
  fid[7,11] <- (1-Fid[7])*(move.primary[7,11]/(move.primary[7,1]+move.primary[7,2]+move.primary[7,3]+move.primary[7,4]+move.primary[7,5]+move.primary[7,6]+move.primary[7,8]+move.primary[7,9]+move.primary[7,10]+move.primary[7,11]+move.primary[7,12]))
  fid[7,12] <- (1-Fid[7])*(move.primary[7,12]/(move.primary[7,1]+move.primary[7,2]+move.primary[7,3]+move.primary[7,4]+move.primary[7,5]+move.primary[7,6]+move.primary[7,8]+move.primary[7,9]+move.primary[7,10]+move.primary[7,11]+move.primary[7,12]))
  fid[7,7] <- Fid[7]
  
  fid[8,1] <- (1-Fid[8])*(move.primary[8,1]/(move.primary[8,1]+move.primary[8,2]+move.primary[8,3]+move.primary[8,4]+move.primary[8,5]+move.primary[8,6]+move.primary[8,7]+move.primary[8,9]+move.primary[8,10]+move.primary[8,11]+move.primary[8,12]))
  fid[8,2] <- (1-Fid[8])*(move.primary[8,2]/(move.primary[8,1]+move.primary[8,2]+move.primary[8,3]+move.primary[8,4]+move.primary[8,5]+move.primary[8,6]+move.primary[8,7]+move.primary[8,9]+move.primary[8,10]+move.primary[8,11]+move.primary[8,12]))
  fid[8,3] <- (1-Fid[8])*(move.primary[8,3]/(move.primary[8,1]+move.primary[8,2]+move.primary[8,3]+move.primary[8,4]+move.primary[8,5]+move.primary[8,6]+move.primary[8,7]+move.primary[8,9]+move.primary[8,10]+move.primary[8,11]+move.primary[8,12]))
  fid[8,4] <- (1-Fid[8])*(move.primary[8,4]/(move.primary[8,1]+move.primary[8,2]+move.primary[8,3]+move.primary[8,4]+move.primary[8,5]+move.primary[8,6]+move.primary[8,7]+move.primary[8,9]+move.primary[8,10]+move.primary[8,11]+move.primary[8,12]))
  fid[8,5] <- (1-Fid[8])*(move.primary[8,5]/(move.primary[8,1]+move.primary[8,2]+move.primary[8,3]+move.primary[8,4]+move.primary[8,5]+move.primary[8,6]+move.primary[8,7]+move.primary[8,9]+move.primary[8,10]+move.primary[8,11]+move.primary[8,12]))
  fid[8,6] <- (1-Fid[8])*(move.primary[8,6]/(move.primary[8,1]+move.primary[8,2]+move.primary[8,3]+move.primary[8,4]+move.primary[8,5]+move.primary[8,6]+move.primary[8,7]+move.primary[8,9]+move.primary[8,10]+move.primary[8,11]+move.primary[8,12]))
  fid[8,7] <- (1-Fid[8])*(move.primary[8,7]/(move.primary[8,1]+move.primary[8,2]+move.primary[8,3]+move.primary[8,4]+move.primary[8,5]+move.primary[8,6]+move.primary[8,7]+move.primary[8,9]+move.primary[8,10]+move.primary[8,11]+move.primary[8,12]))
  fid[8,9] <- (1-Fid[8])*(move.primary[8,9]/(move.primary[8,1]+move.primary[8,2]+move.primary[8,3]+move.primary[8,4]+move.primary[8,5]+move.primary[8,6]+move.primary[8,7]+move.primary[8,9]+move.primary[8,10]+move.primary[8,11]+move.primary[8,12]))
  fid[8,10] <- (1-Fid[8])*(move.primary[8,10]/(move.primary[8,1]+move.primary[8,2]+move.primary[8,3]+move.primary[8,4]+move.primary[8,5]+move.primary[8,6]+move.primary[8,7]+move.primary[8,9]+move.primary[8,10]+move.primary[8,11]+move.primary[8,12]))
  fid[8,11] <- (1-Fid[8])*(move.primary[8,11]/(move.primary[8,1]+move.primary[8,2]+move.primary[8,3]+move.primary[8,4]+move.primary[8,5]+move.primary[8,6]+move.primary[8,7]+move.primary[8,9]+move.primary[8,10]+move.primary[8,11]+move.primary[8,12]))
  fid[8,12] <- (1-Fid[8])*(move.primary[8,12]/(move.primary[8,1]+move.primary[8,2]+move.primary[8,3]+move.primary[8,4]+move.primary[8,5]+move.primary[8,6]+move.primary[8,7]+move.primary[8,9]+move.primary[8,10]+move.primary[8,11]+move.primary[8,12]))
  fid[8,8] <- Fid[8]
  
  fid[9,1] <- (1-Fid[9])*(move.primary[9,1]/(move.primary[9,1]+move.primary[9,2]+move.primary[9,3]+move.primary[9,4]+move.primary[9,5]+move.primary[9,6]+move.primary[9,7]+move.primary[9,8]+move.primary[9,10]+move.primary[9,11]+move.primary[9,12]))
  fid[9,2] <- (1-Fid[9])*(move.primary[9,2]/(move.primary[9,1]+move.primary[9,2]+move.primary[9,3]+move.primary[9,4]+move.primary[9,5]+move.primary[9,6]+move.primary[9,7]+move.primary[9,8]+move.primary[9,10]+move.primary[9,11]+move.primary[9,12]))
  fid[9,3] <- (1-Fid[9])*(move.primary[9,3]/(move.primary[9,1]+move.primary[9,2]+move.primary[9,3]+move.primary[9,4]+move.primary[9,5]+move.primary[9,6]+move.primary[9,7]+move.primary[9,8]+move.primary[9,10]+move.primary[9,11]+move.primary[9,12]))
  fid[9,4] <- (1-Fid[9])*(move.primary[9,4]/(move.primary[9,1]+move.primary[9,2]+move.primary[9,3]+move.primary[9,4]+move.primary[9,5]+move.primary[9,6]+move.primary[9,7]+move.primary[9,8]+move.primary[9,10]+move.primary[9,11]+move.primary[9,12]))
  fid[9,5] <- (1-Fid[9])*(move.primary[9,5]/(move.primary[9,1]+move.primary[9,2]+move.primary[9,3]+move.primary[9,4]+move.primary[9,5]+move.primary[9,6]+move.primary[9,7]+move.primary[9,8]+move.primary[9,10]+move.primary[9,11]+move.primary[9,12]))
  fid[9,6] <- (1-Fid[9])*(move.primary[9,6]/(move.primary[9,1]+move.primary[9,2]+move.primary[9,3]+move.primary[9,4]+move.primary[9,5]+move.primary[9,6]+move.primary[9,7]+move.primary[9,8]+move.primary[9,10]+move.primary[9,11]+move.primary[9,12]))
  fid[9,7] <- (1-Fid[9])*(move.primary[9,7]/(move.primary[9,1]+move.primary[9,2]+move.primary[9,3]+move.primary[9,4]+move.primary[9,5]+move.primary[9,6]+move.primary[9,7]+move.primary[9,8]+move.primary[9,10]+move.primary[9,11]+move.primary[9,12]))
  fid[9,8] <- (1-Fid[9])*(move.primary[9,8]/(move.primary[9,1]+move.primary[9,2]+move.primary[9,3]+move.primary[9,4]+move.primary[9,5]+move.primary[9,6]+move.primary[9,7]+move.primary[9,8]+move.primary[9,10]+move.primary[9,11]+move.primary[9,12]))
  fid[9,10] <- (1-Fid[9])*(move.primary[9,10]/(move.primary[9,1]+move.primary[9,2]+move.primary[9,3]+move.primary[9,4]+move.primary[9,5]+move.primary[9,6]+move.primary[9,7]+move.primary[9,8]+move.primary[9,10]+move.primary[9,11]+move.primary[9,12]))
  fid[9,11] <- (1-Fid[9])*(move.primary[9,11]/(move.primary[9,1]+move.primary[9,2]+move.primary[9,3]+move.primary[9,4]+move.primary[9,5]+move.primary[9,6]+move.primary[9,7]+move.primary[9,8]+move.primary[9,10]+move.primary[9,11]+move.primary[9,12]))
  fid[9,12] <- (1-Fid[9])*(move.primary[9,12]/(move.primary[9,1]+move.primary[9,2]+move.primary[9,3]+move.primary[9,4]+move.primary[9,5]+move.primary[9,6]+move.primary[9,7]+move.primary[9,8]+move.primary[9,10]+move.primary[9,11]+move.primary[9,12]))
  fid[9,9] <- Fid[9]
  
  fid[10,1] <- (1-Fid[10])*(move.primary[10,1]/(move.primary[10,1]+move.primary[10,2]+move.primary[10,3]+move.primary[10,4]+move.primary[10,5]+move.primary[10,6]+move.primary[10,7]+move.primary[10,8]+move.primary[10,9]+move.primary[10,11]+move.primary[10,12]))
  fid[10,2] <- (1-Fid[10])*(move.primary[10,2]/(move.primary[10,1]+move.primary[10,2]+move.primary[10,3]+move.primary[10,4]+move.primary[10,5]+move.primary[10,6]+move.primary[10,7]+move.primary[10,8]+move.primary[10,9]+move.primary[10,11]+move.primary[10,12]))
  fid[10,3] <- (1-Fid[10])*(move.primary[10,3]/(move.primary[10,1]+move.primary[10,2]+move.primary[10,3]+move.primary[10,4]+move.primary[10,5]+move.primary[10,6]+move.primary[10,7]+move.primary[10,8]+move.primary[10,9]+move.primary[10,11]+move.primary[10,12]))
  fid[10,4] <- (1-Fid[10])*(move.primary[10,4]/(move.primary[10,1]+move.primary[10,2]+move.primary[10,3]+move.primary[10,4]+move.primary[10,5]+move.primary[10,6]+move.primary[10,7]+move.primary[10,8]+move.primary[10,9]+move.primary[10,11]+move.primary[10,12]))
  fid[10,5] <- (1-Fid[10])*(move.primary[10,5]/(move.primary[10,1]+move.primary[10,2]+move.primary[10,3]+move.primary[10,4]+move.primary[10,5]+move.primary[10,6]+move.primary[10,7]+move.primary[10,8]+move.primary[10,9]+move.primary[10,11]+move.primary[10,12]))
  fid[10,6] <- (1-Fid[10])*(move.primary[10,6]/(move.primary[10,1]+move.primary[10,2]+move.primary[10,3]+move.primary[10,4]+move.primary[10,5]+move.primary[10,6]+move.primary[10,7]+move.primary[10,8]+move.primary[10,9]+move.primary[10,11]+move.primary[10,12]))
  fid[10,7] <- (1-Fid[10])*(move.primary[10,7]/(move.primary[10,1]+move.primary[10,2]+move.primary[10,3]+move.primary[10,4]+move.primary[10,5]+move.primary[10,6]+move.primary[10,7]+move.primary[10,8]+move.primary[10,9]+move.primary[10,11]+move.primary[10,12]))
  fid[10,8] <- (1-Fid[10])*(move.primary[10,8]/(move.primary[10,1]+move.primary[10,2]+move.primary[10,3]+move.primary[10,4]+move.primary[10,5]+move.primary[10,6]+move.primary[10,7]+move.primary[10,8]+move.primary[10,9]+move.primary[10,11]+move.primary[10,12]))
  fid[10,9] <- (1-Fid[10])*(move.primary[10,9]/(move.primary[10,1]+move.primary[10,2]+move.primary[10,3]+move.primary[10,4]+move.primary[10,5]+move.primary[10,6]+move.primary[10,7]+move.primary[10,8]+move.primary[10,9]+move.primary[10,11]+move.primary[10,12]))
  fid[10,11] <- (1-Fid[10])*(move.primary[10,11]/(move.primary[10,1]+move.primary[10,2]+move.primary[10,3]+move.primary[10,4]+move.primary[10,5]+move.primary[10,6]+move.primary[10,7]+move.primary[10,8]+move.primary[10,9]+move.primary[10,11]+move.primary[10,12]))
  fid[10,12] <- (1-Fid[10])*(move.primary[10,12]/(move.primary[10,1]+move.primary[10,2]+move.primary[10,3]+move.primary[10,4]+move.primary[10,5]+move.primary[10,6]+move.primary[10,7]+move.primary[10,8]+move.primary[10,9]+move.primary[10,11]+move.primary[10,12]))
  fid[10,10] <- Fid[10]
  
  fid[11,1] <- (1-Fid[11])*(move.primary[11,1]/(move.primary[11,1]+move.primary[11,2]+move.primary[11,3]+move.primary[11,4]+move.primary[11,5]+move.primary[11,6]+move.primary[11,7]+move.primary[11,8]+move.primary[11,9]+move.primary[11,10]+move.primary[11,12]))
  fid[11,2] <- (1-Fid[11])*(move.primary[11,2]/(move.primary[11,1]+move.primary[11,2]+move.primary[11,3]+move.primary[11,4]+move.primary[11,5]+move.primary[11,6]+move.primary[11,7]+move.primary[11,8]+move.primary[11,9]+move.primary[11,10]+move.primary[11,12]))
  fid[11,3] <- (1-Fid[11])*(move.primary[11,3]/(move.primary[11,1]+move.primary[11,2]+move.primary[11,3]+move.primary[11,4]+move.primary[11,5]+move.primary[11,6]+move.primary[11,7]+move.primary[11,8]+move.primary[11,9]+move.primary[11,10]+move.primary[11,12]))
  fid[11,4] <- (1-Fid[11])*(move.primary[11,4]/(move.primary[11,1]+move.primary[11,2]+move.primary[11,3]+move.primary[11,4]+move.primary[11,5]+move.primary[11,6]+move.primary[11,7]+move.primary[11,8]+move.primary[11,9]+move.primary[11,10]+move.primary[11,12]))
  fid[11,5] <- (1-Fid[11])*(move.primary[11,5]/(move.primary[11,1]+move.primary[11,2]+move.primary[11,3]+move.primary[11,4]+move.primary[11,5]+move.primary[11,6]+move.primary[11,7]+move.primary[11,8]+move.primary[11,9]+move.primary[11,10]+move.primary[11,12]))
  fid[11,6] <- (1-Fid[11])*(move.primary[11,6]/(move.primary[11,1]+move.primary[11,2]+move.primary[11,3]+move.primary[11,4]+move.primary[11,5]+move.primary[11,6]+move.primary[11,7]+move.primary[11,8]+move.primary[11,9]+move.primary[11,10]+move.primary[11,12]))
  fid[11,7] <- (1-Fid[11])*(move.primary[11,7]/(move.primary[11,1]+move.primary[11,2]+move.primary[11,3]+move.primary[11,4]+move.primary[11,5]+move.primary[11,6]+move.primary[11,7]+move.primary[11,8]+move.primary[11,9]+move.primary[11,10]+move.primary[11,12]))
  fid[11,8] <- (1-Fid[11])*(move.primary[11,8]/(move.primary[11,1]+move.primary[11,2]+move.primary[11,3]+move.primary[11,4]+move.primary[11,5]+move.primary[11,6]+move.primary[11,7]+move.primary[11,8]+move.primary[11,9]+move.primary[11,10]+move.primary[11,12]))
  fid[11,9] <- (1-Fid[11])*(move.primary[11,9]/(move.primary[11,1]+move.primary[11,2]+move.primary[11,3]+move.primary[11,4]+move.primary[11,5]+move.primary[11,6]+move.primary[11,7]+move.primary[11,8]+move.primary[11,9]+move.primary[11,10]+move.primary[11,12]))
  fid[11,10] <- (1-Fid[11])*(move.primary[11,10]/(move.primary[11,1]+move.primary[11,2]+move.primary[11,3]+move.primary[11,4]+move.primary[11,5]+move.primary[11,6]+move.primary[11,7]+move.primary[11,8]+move.primary[11,9]+move.primary[11,10]+move.primary[11,12]))
  fid[11,12] <- (1-Fid[11])*(move.primary[11,12]/(move.primary[11,1]+move.primary[11,2]+move.primary[11,3]+move.primary[11,4]+move.primary[11,5]+move.primary[11,6]+move.primary[11,7]+move.primary[11,8]+move.primary[11,9]+move.primary[11,10]+move.primary[11,12]))
  fid[11,11] <- Fid[11]
  
  fid[12,1] <- (1-Fid[12])*(move.primary[12,1]/(move.primary[12,1]+move.primary[12,2]+move.primary[12,3]+move.primary[12,4]+move.primary[12,5]+move.primary[12,6]+move.primary[12,7]+move.primary[12,8]+move.primary[12,9]+move.primary[12,10]+move.primary[12,11]))
  fid[12,2] <- (1-Fid[12])*(move.primary[12,2]/(move.primary[12,1]+move.primary[12,2]+move.primary[12,3]+move.primary[12,4]+move.primary[12,5]+move.primary[12,6]+move.primary[12,7]+move.primary[12,8]+move.primary[12,9]+move.primary[12,10]+move.primary[12,11]))
  fid[12,3] <- (1-Fid[12])*(move.primary[12,3]/(move.primary[12,1]+move.primary[12,2]+move.primary[12,3]+move.primary[12,4]+move.primary[12,5]+move.primary[12,6]+move.primary[12,7]+move.primary[12,8]+move.primary[12,9]+move.primary[12,10]+move.primary[12,11]))
  fid[12,4] <- (1-Fid[12])*(move.primary[12,4]/(move.primary[12,1]+move.primary[12,2]+move.primary[12,3]+move.primary[12,4]+move.primary[12,5]+move.primary[12,6]+move.primary[12,7]+move.primary[12,8]+move.primary[12,9]+move.primary[12,10]+move.primary[12,11]))
  fid[12,5] <- (1-Fid[12])*(move.primary[12,5]/(move.primary[12,1]+move.primary[12,2]+move.primary[12,3]+move.primary[12,4]+move.primary[12,5]+move.primary[12,6]+move.primary[12,7]+move.primary[12,8]+move.primary[12,9]+move.primary[12,10]+move.primary[12,11]))
  fid[12,6] <- (1-Fid[12])*(move.primary[12,6]/(move.primary[12,1]+move.primary[12,2]+move.primary[12,3]+move.primary[12,4]+move.primary[12,5]+move.primary[12,6]+move.primary[12,7]+move.primary[12,8]+move.primary[12,9]+move.primary[12,10]+move.primary[12,11]))
  fid[12,7] <- (1-Fid[12])*(move.primary[12,7]/(move.primary[12,1]+move.primary[12,2]+move.primary[12,3]+move.primary[12,4]+move.primary[12,5]+move.primary[12,6]+move.primary[12,7]+move.primary[12,8]+move.primary[12,9]+move.primary[12,10]+move.primary[12,11]))
  fid[12,8] <- (1-Fid[12])*(move.primary[12,8]/(move.primary[12,1]+move.primary[12,2]+move.primary[12,3]+move.primary[12,4]+move.primary[12,5]+move.primary[12,6]+move.primary[12,7]+move.primary[12,8]+move.primary[12,9]+move.primary[12,10]+move.primary[12,11]))
  fid[12,9] <- (1-Fid[12])*(move.primary[12,9]/(move.primary[12,1]+move.primary[12,2]+move.primary[12,3]+move.primary[12,4]+move.primary[12,5]+move.primary[12,6]+move.primary[12,7]+move.primary[12,8]+move.primary[12,9]+move.primary[12,10]+move.primary[12,11]))
  fid[12,10] <- (1-Fid[12])*(move.primary[12,10]/(move.primary[12,1]+move.primary[12,2]+move.primary[12,3]+move.primary[12,4]+move.primary[12,5]+move.primary[12,6]+move.primary[12,7]+move.primary[12,8]+move.primary[12,9]+move.primary[12,10]+move.primary[12,11]))
  fid[12,11] <- (1-Fid[12])*(move.primary[12,11]/(move.primary[12,1]+move.primary[12,2]+move.primary[12,3]+move.primary[12,4]+move.primary[12,5]+move.primary[12,6]+move.primary[12,7]+move.primary[12,8]+move.primary[12,9]+move.primary[12,10]+move.primary[12,11]))
  fid[12,12] <- Fid[12]
  
  # transition probability matrices between periods
  gamma <- array(0, dim=c(14, 14, T-1))
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
  
  gamma[2,2,] <- s * fid[1,1] # Must survive and then move (or not move) between years
  gamma[2,3,] <- s * fid[1,2]
  gamma[2,4,] <- s * fid[1,3]
  gamma[2,5,] <- s * fid[1,4]
  gamma[2,6,] <- s * fid[1,5]
  gamma[2,7,] <- s * fid[1,6]
  gamma[2,8,] <- s * fid[1,7]
  gamma[2,9,] <- s * fid[1,8]
  gamma[2,10,] <- s * fid[1,9]
  gamma[2,11,] <- s * fid[1,10]
  gamma[2,12,] <- s * fid[1,11]
  gamma[2,13,] <- s * fid[1,12]
  gamma[2,14,] <- 1-s
  
  gamma[3,2,] <- s * fid[2,1]
  gamma[3,3,] <- s * fid[2,2]
  gamma[3,4,] <- s * fid[2,3]
  gamma[3,5,] <- s * fid[2,4]
  gamma[3,6,] <- s * fid[2,5]
  gamma[3,7,] <- s * fid[2,6]
  gamma[3,8,] <- s * fid[2,7]
  gamma[3,9,] <- s * fid[2,8]
  gamma[3,10,] <- s * fid[2,9]
  gamma[3,11,] <- s * fid[2,10]
  gamma[3,12,] <- s * fid[2,11]
  gamma[3,13,] <- s * fid[2,12]
  gamma[3,14,] <- 1-s
  
  gamma[4,2,] <- s * fid[3,1]
  gamma[4,3,] <- s * fid[3,2]
  gamma[4,4,] <- s * fid[3,3]
  gamma[4,5,] <- s * fid[3,4]
  gamma[4,6,] <- s * fid[3,5]
  gamma[4,7,] <- s * fid[3,6]
  gamma[4,8,] <- s * fid[3,7]
  gamma[4,9,] <- s * fid[3,8]
  gamma[4,10,] <- s * fid[3,9]
  gamma[4,11,] <- s * fid[3,10]
  gamma[4,12,] <- s * fid[3,11]
  gamma[4,13,] <- s * fid[3,12]
  gamma[4,14,] <- 1-s
  
  gamma[5,2,] <- s * fid[4,1]
  gamma[5,3,] <- s * fid[4,2]
  gamma[5,4,] <- s * fid[4,3]
  gamma[5,5,] <- s * fid[4,4]
  gamma[5,6,] <- s * fid[4,5]
  gamma[5,7,] <- s * fid[4,6]
  gamma[5,8,] <- s * fid[4,7]
  gamma[5,9,] <- s * fid[4,8]
  gamma[5,10,] <- s * fid[4,9]
  gamma[5,11,] <- s * fid[4,10]
  gamma[5,12,] <- s * fid[4,11]
  gamma[5,13,] <- s * fid[4,12]
  gamma[5,14,] <- 1-s
  
  gamma[6,2,] <- s * fid[5,1]
  gamma[6,3,] <- s * fid[5,2]
  gamma[6,4,] <- s * fid[5,3]
  gamma[6,5,] <- s * fid[5,4]
  gamma[6,6,] <- s * fid[5,5]
  gamma[6,7,] <- s * fid[5,6]
  gamma[6,8,] <- s * fid[5,7]
  gamma[6,9,] <- s * fid[5,8]
  gamma[6,10,] <- s * fid[5,9]
  gamma[6,11,] <- s * fid[5,10]
  gamma[6,12,] <- s * fid[5,11]
  gamma[6,13,] <- s * fid[5,12]
  gamma[6,14,] <- 1-s
  
  gamma[7,2,] <- s * fid[6,1]
  gamma[7,3,] <- s * fid[6,2]
  gamma[7,4,] <- s * fid[6,3]
  gamma[7,5,] <- s * fid[6,4]
  gamma[7,6,] <- s * fid[6,5]
  gamma[7,7,] <- s * fid[6,6]
  gamma[7,8,] <- s * fid[6,7]
  gamma[7,9,] <- s * fid[6,8]
  gamma[7,10,] <- s * fid[6,9]
  gamma[7,11,] <- s * fid[6,10]
  gamma[7,12,] <- s * fid[6,11]
  gamma[7,13,] <- s * fid[6,12]
  gamma[7,14,] <- 1-s
  
  gamma[8,2,] <- s * fid[7,1]
  gamma[8,3,] <- s * fid[7,2]
  gamma[8,4,] <- s * fid[7,3]
  gamma[8,5,] <- s * fid[7,4]
  gamma[8,6,] <- s * fid[7,5]
  gamma[8,7,] <- s * fid[7,6]
  gamma[8,8,] <- s * fid[7,7]
  gamma[8,9,] <- s * fid[7,8]
  gamma[8,10,] <- s * fid[7,9]
  gamma[8,11,] <- s * fid[7,10]
  gamma[8,12,] <- s * fid[7,11]
  gamma[8,13,] <- s * fid[7,12]
  gamma[8,14,] <- 1-s
  
  gamma[9,2,] <- s * fid[8,1]
  gamma[9,3,] <- s * fid[8,2]
  gamma[9,4,] <- s * fid[8,3]
  gamma[9,5,] <- s * fid[8,4]
  gamma[9,6,] <- s * fid[8,5]
  gamma[9,7,] <- s * fid[8,6]
  gamma[9,8,] <- s * fid[8,7]
  gamma[9,9,] <- s * fid[8,8]
  gamma[9,10,] <- s * fid[8,9]
  gamma[9,11,] <- s * fid[8,10]
  gamma[9,12,] <- s * fid[8,11]
  gamma[9,13,] <- s * fid[8,12]
  gamma[9,14,] <- 1-s
  
  gamma[10,2,] <- s * fid[9,1]
  gamma[10,3,] <- s * fid[9,2]
  gamma[10,4,] <- s * fid[9,3]
  gamma[10,5,] <- s * fid[9,4]
  gamma[10,6,] <- s * fid[9,5]
  gamma[10,7,] <- s * fid[9,6]
  gamma[10,8,] <- s * fid[9,7]
  gamma[10,9,] <- s * fid[9,8]
  gamma[10,10,] <- s * fid[9,9]
  gamma[10,11,] <- s * fid[9,10]
  gamma[10,12,] <- s * fid[9,11]
  gamma[10,13,] <- s * fid[9,12]
  gamma[10,14,] <- 1-s
  
  gamma[11,2,] <- s * fid[10,1]
  gamma[11,3,] <- s * fid[10,2]
  gamma[11,4,] <- s * fid[10,3]
  gamma[11,5,] <- s * fid[10,4]
  gamma[11,6,] <- s * fid[10,5]
  gamma[11,7,] <- s * fid[10,6]
  gamma[11,8,] <- s * fid[10,7]
  gamma[11,9,] <- s * fid[10,8]
  gamma[11,10,] <- s * fid[10,9]
  gamma[11,11,] <- s * fid[10,10]
  gamma[11,12,] <- s * fid[10,11]
  gamma[11,13,] <- s * fid[10,12]
  gamma[11,14,] <- 1-s
  
  gamma[12,2,] <- s * fid[11,1]
  gamma[12,3,] <- s * fid[11,2]
  gamma[12,4,] <- s * fid[11,3]
  gamma[12,5,] <- s * fid[11,4]
  gamma[12,6,] <- s * fid[11,5]
  gamma[12,7,] <- s * fid[11,6]
  gamma[12,8,] <- s * fid[11,7]
  gamma[12,9,] <- s * fid[11,8]
  gamma[12,10,] <- s * fid[11,9]
  gamma[12,11,] <- s * fid[11,10]
  gamma[12,12,] <- s * fid[11,11]
  gamma[12,13,] <- s * fid[11,12]
  gamma[12,14,] <- 1-s
  
  gamma[13,2,] <- s * fid[12,1]
  gamma[13,3,] <- s * fid[12,2]
  gamma[13,4,] <- s * fid[12,3]
  gamma[13,5,] <- s * fid[12,4]
  gamma[13,6,] <- s * fid[12,5]
  gamma[13,7,] <- s * fid[12,6]
  gamma[13,8,] <- s * fid[12,7]
  gamma[13,9,] <- s * fid[12,8]
  gamma[13,10,] <- s * fid[12,9]
  gamma[13,11,] <- s * fid[12,10]
  gamma[13,12,] <- s * fid[12,11]
  gamma[13,13,] <- s * fid[12,12]
  gamma[13,14,] <- 1-s
  
  gamma[14,14,] <- 1 # stay dead
  
  
  
  # transition probability matrices within periods
  gammat <- NULL
  for (t in 1:T)  {
    gammat[[t]] <- array(0, dim=c(14, 14, K[t]-1))
    gammat[[t]][1,1,] <- 1 - betastar[[t]][2:K[t]]
    gammat[[t]][1,2,] <- betastar[[t]][2:K[t]]*alpha[1]
    gammat[[t]][1,3,] <- betastar[[t]][2:K[t]]*alpha[2]
    gammat[[t]][1,4,] <- betastar[[t]][2:K[t]]*alpha[3]
    gammat[[t]][1,5,] <- betastar[[t]][2:K[t]]*alpha[4]
    gammat[[t]][1,6,] <- betastar[[t]][2:K[t]]*alpha[5]
    gammat[[t]][1,7,] <- betastar[[t]][2:K[t]]*alpha[6]
    gammat[[t]][1,8,] <- betastar[[t]][2:K[t]]*alpha[7]
    gammat[[t]][1,9,] <- betastar[[t]][2:K[t]]*alpha[8]
    gammat[[t]][1,10,] <- betastar[[t]][2:K[t]]*alpha[9]
    gammat[[t]][1,11,] <- betastar[[t]][2:K[t]]*alpha[10]
    gammat[[t]][1,12,] <- betastar[[t]][2:K[t]]*alpha[11]
    gammat[[t]][1,13,] <- betastar[[t]][2:K[t]]*alpha[12]
    
    gammat[[t]][2,2,] <- phistar[[t]]*psi[1,1]
    gammat[[t]][2,3,] <- phistar[[t]]*psi[1,2]
    gammat[[t]][2,4,] <- phistar[[t]]*psi[1,3]
    gammat[[t]][2,5,] <- phistar[[t]]*psi[1,4]
    gammat[[t]][2,6,] <- phistar[[t]]*psi[1,5]
    gammat[[t]][2,7,] <- phistar[[t]]*psi[1,6]
    gammat[[t]][2,8,] <- phistar[[t]]*psi[1,7]
    gammat[[t]][2,9,] <- phistar[[t]]*psi[1,8]
    gammat[[t]][2,10,] <- phistar[[t]]*psi[1,9]
    gammat[[t]][2,11,] <- phistar[[t]]*psi[1,10]
    gammat[[t]][2,12,] <- phistar[[t]]*psi[1,11]
    gammat[[t]][2,13,] <- phistar[[t]]*psi[1,12]
    gammat[[t]][2,14,] <- 1-phistar[[t]]
    
    gammat[[t]][3,2,] <- phistar[[t]]*psi[2,1]
    gammat[[t]][3,3,] <- phistar[[t]]*psi[2,2]
    gammat[[t]][3,4,] <- phistar[[t]]*psi[2,3]
    gammat[[t]][3,5,] <- phistar[[t]]*psi[2,4]
    gammat[[t]][3,6,] <- phistar[[t]]*psi[2,5]
    gammat[[t]][3,7,] <- phistar[[t]]*psi[2,6]
    gammat[[t]][3,8,] <- phistar[[t]]*psi[2,7]
    gammat[[t]][3,9,] <- phistar[[t]]*psi[2,8]
    gammat[[t]][3,10,] <- phistar[[t]]*psi[2,9]
    gammat[[t]][3,11,] <- phistar[[t]]*psi[2,10]
    gammat[[t]][3,12,] <- phistar[[t]]*psi[2,11]
    gammat[[t]][3,13,] <- phistar[[t]]*psi[2,12]
    gammat[[t]][3,14,] <- 1-phistar[[t]]
    
    gammat[[t]][4,2,] <- phistar[[t]]*psi[3,1]
    gammat[[t]][4,3,] <- phistar[[t]]*psi[3,2]
    gammat[[t]][4,4,] <- phistar[[t]]*psi[3,3]
    gammat[[t]][4,5,] <- phistar[[t]]*psi[3,4]
    gammat[[t]][4,6,] <- phistar[[t]]*psi[3,5]
    gammat[[t]][4,7,] <- phistar[[t]]*psi[3,6]
    gammat[[t]][4,8,] <- phistar[[t]]*psi[3,7]
    gammat[[t]][4,9,] <- phistar[[t]]*psi[3,8]
    gammat[[t]][4,10,] <- phistar[[t]]*psi[3,9]
    gammat[[t]][4,11,] <- phistar[[t]]*psi[3,10]
    gammat[[t]][4,12,] <- phistar[[t]]*psi[3,11]
    gammat[[t]][4,13,] <- phistar[[t]]*psi[3,12]
    gammat[[t]][4,14,] <- 1-phistar[[t]]
    
    gammat[[t]][5,2,] <- phistar[[t]]*psi[4,1]
    gammat[[t]][5,3,] <- phistar[[t]]*psi[4,2]
    gammat[[t]][5,4,] <- phistar[[t]]*psi[4,3]
    gammat[[t]][5,5,] <- phistar[[t]]*psi[4,4]
    gammat[[t]][5,6,] <- phistar[[t]]*psi[4,5]
    gammat[[t]][5,7,] <- phistar[[t]]*psi[4,6]
    gammat[[t]][5,8,] <- phistar[[t]]*psi[4,7]
    gammat[[t]][5,9,] <- phistar[[t]]*psi[4,8]
    gammat[[t]][5,10,] <- phistar[[t]]*psi[4,9]
    gammat[[t]][5,11,] <- phistar[[t]]*psi[4,10]
    gammat[[t]][5,12,] <- phistar[[t]]*psi[4,11]
    gammat[[t]][5,13,] <- phistar[[t]]*psi[4,12]
    gammat[[t]][5,14,] <- 1-phistar[[t]]
    
    gammat[[t]][6,2,] <- phistar[[t]]*psi[5,1]
    gammat[[t]][6,3,] <- phistar[[t]]*psi[5,2]
    gammat[[t]][6,4,] <- phistar[[t]]*psi[5,3]
    gammat[[t]][6,5,] <- phistar[[t]]*psi[5,4]
    gammat[[t]][6,6,] <- phistar[[t]]*psi[5,5]
    gammat[[t]][6,7,] <- phistar[[t]]*psi[5,6]
    gammat[[t]][6,8,] <- phistar[[t]]*psi[5,7]
    gammat[[t]][6,9,] <- phistar[[t]]*psi[5,8]
    gammat[[t]][6,10,] <- phistar[[t]]*psi[5,9]
    gammat[[t]][6,11,] <- phistar[[t]]*psi[5,10]
    gammat[[t]][6,12,] <- phistar[[t]]*psi[5,11]
    gammat[[t]][6,13,] <- phistar[[t]]*psi[5,12]
    gammat[[t]][6,14,] <- 1-phistar[[t]]
    
    gammat[[t]][7,2,] <- phistar[[t]]*psi[6,1]
    gammat[[t]][7,3,] <- phistar[[t]]*psi[6,2]
    gammat[[t]][7,4,] <- phistar[[t]]*psi[6,3]
    gammat[[t]][7,5,] <- phistar[[t]]*psi[6,4]
    gammat[[t]][7,6,] <- phistar[[t]]*psi[6,5]
    gammat[[t]][7,7,] <- phistar[[t]]*psi[6,6]
    gammat[[t]][7,8,] <- phistar[[t]]*psi[6,7]
    gammat[[t]][7,9,] <- phistar[[t]]*psi[6,8]
    gammat[[t]][7,10,] <- phistar[[t]]*psi[6,9]
    gammat[[t]][7,11,] <- phistar[[t]]*psi[6,10]
    gammat[[t]][7,12,] <- phistar[[t]]*psi[6,11]
    gammat[[t]][7,13,] <- phistar[[t]]*psi[6,12]
    gammat[[t]][7,14,] <- 1-phistar[[t]]
    
    gammat[[t]][8,2,] <- phistar[[t]]*psi[7,1]
    gammat[[t]][8,3,] <- phistar[[t]]*psi[7,2]
    gammat[[t]][8,4,] <- phistar[[t]]*psi[7,3]
    gammat[[t]][8,5,] <- phistar[[t]]*psi[7,4]
    gammat[[t]][8,6,] <- phistar[[t]]*psi[7,5]
    gammat[[t]][8,7,] <- phistar[[t]]*psi[7,6]
    gammat[[t]][8,8,] <- phistar[[t]]*psi[7,7]
    gammat[[t]][8,9,] <- phistar[[t]]*psi[7,8]
    gammat[[t]][8,10,] <- phistar[[t]]*psi[7,9]
    gammat[[t]][8,11,] <- phistar[[t]]*psi[7,10]
    gammat[[t]][8,12,] <- phistar[[t]]*psi[7,11]
    gammat[[t]][8,13,] <- phistar[[t]]*psi[7,12]
    gammat[[t]][8,14,] <- 1-phistar[[t]]
    
    gammat[[t]][9,2,] <- phistar[[t]]*psi[8,1]
    gammat[[t]][9,3,] <- phistar[[t]]*psi[8,2]
    gammat[[t]][9,4,] <- phistar[[t]]*psi[8,3]
    gammat[[t]][9,5,] <- phistar[[t]]*psi[8,4]
    gammat[[t]][9,6,] <- phistar[[t]]*psi[8,5]
    gammat[[t]][9,7,] <- phistar[[t]]*psi[8,6]
    gammat[[t]][9,8,] <- phistar[[t]]*psi[8,7]
    gammat[[t]][9,9,] <- phistar[[t]]*psi[8,8]
    gammat[[t]][9,10,] <- phistar[[t]]*psi[8,9]
    gammat[[t]][9,11,] <- phistar[[t]]*psi[8,10]
    gammat[[t]][9,12,] <- phistar[[t]]*psi[8,11]
    gammat[[t]][9,13,] <- phistar[[t]]*psi[8,12]
    gammat[[t]][9,14,] <- 1-phistar[[t]]
    
    gammat[[t]][10,2,] <- phistar[[t]]*psi[9,1]
    gammat[[t]][10,3,] <- phistar[[t]]*psi[9,2]
    gammat[[t]][10,4,] <- phistar[[t]]*psi[9,3]
    gammat[[t]][10,5,] <- phistar[[t]]*psi[9,4]
    gammat[[t]][10,6,] <- phistar[[t]]*psi[9,5]
    gammat[[t]][10,7,] <- phistar[[t]]*psi[9,6]
    gammat[[t]][10,8,] <- phistar[[t]]*psi[9,7]
    gammat[[t]][10,9,] <- phistar[[t]]*psi[9,8]
    gammat[[t]][10,10,] <- phistar[[t]]*psi[9,9]
    gammat[[t]][10,11,] <- phistar[[t]]*psi[9,10]
    gammat[[t]][10,12,] <- phistar[[t]]*psi[9,11]
    gammat[[t]][10,13,] <- phistar[[t]]*psi[9,12]
    gammat[[t]][10,14,] <- 1-phistar[[t]]
    
    gammat[[t]][11,2,] <- phistar[[t]]*psi[10,1]
    gammat[[t]][11,3,] <- phistar[[t]]*psi[10,2]
    gammat[[t]][11,4,] <- phistar[[t]]*psi[10,3]
    gammat[[t]][11,5,] <- phistar[[t]]*psi[10,4]
    gammat[[t]][11,6,] <- phistar[[t]]*psi[10,5]
    gammat[[t]][11,7,] <- phistar[[t]]*psi[10,6]
    gammat[[t]][11,8,] <- phistar[[t]]*psi[10,7]
    gammat[[t]][11,9,] <- phistar[[t]]*psi[10,8]
    gammat[[t]][11,10,] <- phistar[[t]]*psi[10,9]
    gammat[[t]][11,11,] <- phistar[[t]]*psi[10,10]
    gammat[[t]][11,12,] <- phistar[[t]]*psi[10,11]
    gammat[[t]][11,13,] <- phistar[[t]]*psi[10,12]
    gammat[[t]][11,14,] <- 1-phistar[[t]]
    
    gammat[[t]][12,2,] <- phistar[[t]]*psi[11,1]
    gammat[[t]][12,3,] <- phistar[[t]]*psi[11,2]
    gammat[[t]][12,4,] <- phistar[[t]]*psi[11,3]
    gammat[[t]][12,5,] <- phistar[[t]]*psi[11,4]
    gammat[[t]][12,6,] <- phistar[[t]]*psi[11,5]
    gammat[[t]][12,7,] <- phistar[[t]]*psi[11,6]
    gammat[[t]][12,8,] <- phistar[[t]]*psi[11,7]
    gammat[[t]][12,9,] <- phistar[[t]]*psi[11,8]
    gammat[[t]][12,10,] <- phistar[[t]]*psi[11,9]
    gammat[[t]][12,11,] <- phistar[[t]]*psi[11,10]
    gammat[[t]][12,12,] <- phistar[[t]]*psi[11,11]
    gammat[[t]][12,13,] <- phistar[[t]]*psi[11,12]
    gammat[[t]][12,14,] <- 1-phistar[[t]]
    
    gammat[[t]][13,2,] <- phistar[[t]]*psi[12,1]
    gammat[[t]][13,3,] <- phistar[[t]]*psi[12,2]
    gammat[[t]][13,4,] <- phistar[[t]]*psi[12,3]
    gammat[[t]][13,5,] <- phistar[[t]]*psi[12,4]
    gammat[[t]][13,6,] <- phistar[[t]]*psi[12,5]
    gammat[[t]][13,7,] <- phistar[[t]]*psi[12,6]
    gammat[[t]][13,8,] <- phistar[[t]]*psi[12,7]
    gammat[[t]][13,9,] <- phistar[[t]]*psi[12,8]
    gammat[[t]][13,10,] <- phistar[[t]]*psi[12,9]
    gammat[[t]][13,11,] <- phistar[[t]]*psi[12,10]
    gammat[[t]][13,12,] <- phistar[[t]]*psi[12,11]
    gammat[[t]][13,13,] <- phistar[[t]]*psi[12,12]
    gammat[[t]][13,14,] <- 1-phistar[[t]]
    
    gammat[[t]][14,14,] <- 1
  }
  
  # observation matrices within periods
  pmat <- NULL
  for (t in 1:T)  {
    pmat[[t]] <- array(0, dim=c(14,14,13,K[t]))
    pmat[[t]][1,1,1,] <- 1
    pmat[[t]][2,2,1,] <- 1-p[[t]]
    pmat[[t]][3,3,1,] <- 1-p[[t]]
    pmat[[t]][4,4,1,] <- 1-p[[t]]
    pmat[[t]][5,5,1,] <- 1-p[[t]]
    pmat[[t]][6,6,1,] <- 1-p[[t]]
    pmat[[t]][7,7,1,] <- 1-p[[t]]
    pmat[[t]][8,8,1,] <- 1-p[[t]]
    pmat[[t]][9,9,1,] <- 1-p[[t]]
    pmat[[t]][10,10,1,] <- 1-p[[t]]
    pmat[[t]][11,11,1,] <- 1-p[[t]]
    pmat[[t]][12,12,1,] <- 1-p[[t]]
    pmat[[t]][13,13,1,] <- 1-p[[t]]
    pmat[[t]][14,14,1,] <- 1
    
    pmat[[t]][2,2,2,] <- p[[t]]
    pmat[[t]][3,3,3,] <- p[[t]]
    pmat[[t]][4,4,4,] <- p[[t]]
    pmat[[t]][5,5,5,] <- p[[t]]
    pmat[[t]][6,6,6,] <- p[[t]]
    pmat[[t]][7,7,7,] <- p[[t]]
    pmat[[t]][8,8,8,] <- p[[t]]
    pmat[[t]][9,9,9,] <- p[[t]]
    pmat[[t]][10,10,10,] <- p[[t]]
    pmat[[t]][11,11,11,] <- p[[t]]
    pmat[[t]][12,12,12,] <- p[[t]]
    pmat[[t]][13,13,13,] <- p[[t]]
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
  
  # state probabilities for occasion 1 within periods
  pionet <- NULL
  for (t in 1:T)  {
    pionet[[t]] <- rep(0,14)
    pionet[[t]][1] <- 1-betastar[[t]][1]
    pionet[[t]][2] <- betastar[[t]][1]*alpha[1]
    pionet[[t]][3] <- betastar[[t]][1]*alpha[2]
    pionet[[t]][4] <- betastar[[t]][1]*alpha[3]
    pionet[[t]][5] <- betastar[[t]][1]*alpha[4]
    pionet[[t]][6] <- betastar[[t]][1]*alpha[5]
    pionet[[t]][7] <- betastar[[t]][1]*alpha[6]
    pionet[[t]][8] <- betastar[[t]][1]*alpha[7]
    pionet[[t]][9] <- betastar[[t]][1]*alpha[8]
    pionet[[t]][10] <- betastar[[t]][1]*alpha[9]
    pionet[[t]][11] <- betastar[[t]][1]*alpha[10]
    pionet[[t]][12] <- betastar[[t]][1]*alpha[11]
    pionet[[t]][13] <- betastar[[t]][1]*alpha[12]
  }
  
  # return all parameters
  return(list('n.missed'=n.missed, 'N'=N, 'r'=r, 'rstar'=rstar, 's'=s, 'betaintercept'=betaintercept, 'betagradient'=betagradient, 'betalogit'=betalogit, 'betastar'=betastar, 'alpha'=alpha, 'p'=p, 'sigma' = sigma, 'sigma.primary' = sigma.primary, 'psi'=psi, 'fid' = fid, 'phiintercept'= phiintercept, 'phigradient'=phigradient, 'philogit'=philogit, 'phistar'=phistar, 'gamma'=gamma, 'gammat'=gammat, 'pmat'=pmat, 'pione'=pione, 'pionet'=pionet))
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
    test <- test%*%matrix(1, nrow=14, ncol=1)
    likzerot[t] <- test
  }
  
  # observation matrix for missed individuals
  pzero <- array(0, dim=c(14,14,T))
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
  
  # unobserved individuals 
  likzero <- params$pione
  likzero <- likzero%*%pzero[,,1]
  for (t in 1:(T-1))  {
    likzero <- likzero%*%params$gamma[,,t]
    likzero <- likzero%*%pzero[,,t+1]
  }
  likzero <- likzero%*%matrix(1, nrow=14, ncol=1)
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
        likobst[t] <- test%*%matrix(1, nrow=14, ncol=1)
      }
    }
    
    # observation matrices
    pobs <- array(0, dim=c(14,14,13,T))
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
    likobs[i] <- test%*%matrix(1, nrow=14, ncol=1)
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
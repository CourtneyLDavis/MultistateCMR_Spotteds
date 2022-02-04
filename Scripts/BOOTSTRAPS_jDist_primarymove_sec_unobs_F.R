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
# This script defines the function to bootstrap the capture histories, 
# perform the optimization of the HMM likelihood and store the parameter 
# estimates at each bootstrap iteration.
# -------------------------------------------------------------------------------

# -------------------------------------------------------------------------------
# Inputs: param - vector of parameters to be passed to param.unpack
#         n - number of observed individuals
#         T - number of years
#         K - number of occasions in each year
#         X - list of capture histories
#         Z - attendance history
#         n.hist - frequency of histories
#         distSq - matrix describing distances between sites
#         nboot - the number of bootstraps to perform
# Outputs: matrix of parameters estimates
#          !!! WRITES FILES OUT AFTER EACH BOOTSTRAP !!!      
# -------------------------------------------------------------------------------

bootstrap <- function(param, n, T, K, X, Z, n.hist, distSq, nboot)  {
  
  # indicator of number of successful bootstraps
  flag <- 0
  
  # storage matrix
  res <- matrix(0, nrow=nboot, ncol=length(param))
  
  # names for files
  Znames <- rep(0, nboot)
  Xnames <- matrix(0, nrow=nboot, ncol=T)
  for (i in 1:nboot)  {
    Znames[i] <- paste('Zboot', i, '.txt', sep='')
    for (t in 1:T)  {
      Xnames[i,t] <- paste('Xboot', i, 'year', t, '.txt', sep='')
    }
  }
  
  # perform the bootstrap
  while (flag < nboot)  {
    
    # bootstrap data by sampling individuals
    bootsample <- sample(1:n, n, replace=T)
    Z.boot <- Z[bootsample,]
    X.boot <- NULL
    for (t in 1:T)  {
      X.boot[[t]] <- X[[t]][bootsample,]
    }
    
    write.table(Z.boot, Znames[flag+1], row.names=F, col.names=F)
    for (t in 1:T)  {
      write.table(X.boot[[t]], Xnames[(flag+1),t], row.names=F, col.names=F)
    }
    
    # run optimiser to maximise HMM likelihood
    opt <- nlm(likelihood, param, n, T, K, X.boot, Z.boot, n.hist, distSq)
    if (opt$code == 1 | opt$code == 2)  {
      flag2 <- 3
    } else  {
      flag2 <- 1
    }
    
    # check for convergence
    while (flag2 != 3)  {
      opt <- nlm(likelihood, opt$estimate, n, T, K, X.boot, Z.boot, n.hist, distSq)
      if (opt$code == 1 | opt$code == 2)  {
        flag2 <- 3
      } else  {
        flag2 <- flag2+1
      }
    }
    
    # store estimates
    res[flag+1,] <- opt$estimate
    
    # increase bootstrap number
    if (flag2 == 3)  {
      flag <- flag+1
    }
    
    # print and write to file
    print(c('bootstrap', flag, 'finished'))
    write.table(res, 'bootstrapres.txt', row.names=F, col.names=F)
  }
  
  # results
  return(res)
}
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

source('Scripts/FUNCTIONS_jDist_fullmodel_M.R')

# read in data for primary level (between-year observations)
Z <- read.table('Data/attendance_12ponds_primarymove_males.txt', header=F)

# read in data for secondary level (within-year observations)
X <- NULL
X[[1]] <- read.table('Data/state2014_12ponds_males.txt', header=F)
X[[2]] <- read.table('Data/state2015_12ponds_males.txt', header=F)
X[[3]] <- read.table('Data/state2016_12ponds_males.txt', header=F)
X[[4]] <- read.table('Data/state2017_12ponds_males.txt', header=F)
X[[5]] <- read.table('Data/state2018_12ponds_males.txt', header=F)
X[[6]] <- read.table('Data/state2019_12ponds_males.txt', header=F)

# read in distance matrix and convert to kilometers
distSq <- read.csv("Data/DistanceMatrix_AMMACMR.csv", header = F)
colnames(distSq) <- NULL
distSq <- as.matrix(distSq)
distSq <- distSq/1000

# set constants
n <- length(X[[1]][,1])
T <- 6
K <- rep(0, T)
for (t in 1:T)  {
  K[t] <- length(X[[t]][1,])
}
n.hist <- rep(1,n) # data contain no repeat histories

# set initial values for all parameters
param <- c(2.0914, 
           0.11346734, -1.386294, -0.4054651, -0.8472979, -0.67278778, 
           1.47823538, 
           0.92719425, 
           -0.93515200, -1.72337910, -0.91372267, -4.04938570, -1.26881758, -2.67502349,
           rep(-1.098012, 11),
           rep(1.47823538,6),
           -2,
           rep(-1.09, 12),
           -1.40,
           -1.50819940, -10.34739280, -2.04477353, -1.72637293, -13.12783719, -1.47074036,
           1.234, 
           rep(-1.09,12))

# check the vector is unpacked correctly
params <- param.unpack(param, n, T, K, distSq)

# try likelihood function at starting values
lik <- likelihood(param, n, T, K, X, Z, n.hist, distSq)

# optimise likelihood
opt <- nlm(likelihood, param, n, T, K, X, Z, n.hist, distSq, iterlim = 500)

# MLEs (using output from nlm)
res <- param.unpack(opt$estimate, n, T, K, distSq)


# run to calculate bootstrap confidence intervals
source('BOOTSTRAPS_jDist_fullmodel_M.R')

# run bootstrap (using output from nlm)
resboot <- bootstrap(opt$estimate, n, T, K, X, Z, n.hist, distSq, 50)

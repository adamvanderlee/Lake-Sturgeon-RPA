################################################################################
#        1         2         3         4         5          6        7         8                   
#2345678901234567890123456789012345678901234567890123456789012345678901234567890
################################################################################
# Matrix functions
#-------------------------------------------------------------------------------

## For substituting variables into an expression of variables
## functions courtesy of Gabor Grothendieck from R help forum
esub <- function(expr, sublist) do.call("substitute", list(expr, sublist))

proc <- function(e, env = parent.frame()) {
   for(nm in all.vars(e)) {
      if (exists(nm, env) && is.language(g <- get(nm, env))) {
         if (is.expression(g)) g <- g[[1]]
            g <- Recall(g, env)
            L <- list(g)
            names(L) <- nm
             e <- esub(e, L)
      }
    }
     e
}
#

#-------------------------------------------------------------------------------
# For parameterizing projection matrix from expression 

pmx_eval = function(mx, vars){
  matrix(sapply(mx, eval, vars),  sqrt(length(mx)), sqrt(length(mx)), byrow = T)
}

#-------------------------------------------------------------------------------
# Expression functions
# Create a vector of age-specific mean values from a series of expressions
# Evalute mean values of vital rates from expressions that share a suffix with
# Differeing prefixes (i.e. s0, S1, etc.)

# Creates vector of names expression that share a prefix
Ex_vec <- function(prefix, suffix){
  Exp <- sapply(suffix, function(x) get(paste0(prefix, x))) # Create vecotr of expression
  names(Exp) <- paste0(prefix, suffix)                      # Name expression
  Exp
}

# substitute values into expression and evaluate
Ex_means <- function(prefix, suffix, data) {
  sapply(as.expression(sapply(Ex_vec(prefix, suffix), proc)), eval, data)
}

#-------------------------------------------------------------------------------
# vec-permutation matrix
# function to derive the vec-permutation matrix used to combine the demographic
# projection matrix and dispersal matrix of a metapopulation when the dpm is
# organized by patches and the dm is organized by stages (see Hunter and Caswell 2005)
# output is maxtrix P with dimesions sp X sp.
# vec-permutation matrix is calc as the kronecker product of E and the transpose of E
# sumed over i stages and j patches, where E is a s by p sized matrix with
# a 1 at element E[i,j]
VPM <- function(s, # number of stages in projection matrix
                p) # number of patches in projection matrix
{
    
    # initialize vec-permutation matrix dimesions sp X sp
    P = matrix(0, nrow = s*p, ncol = s*p) 
    
    # Est vec-permutation matrix
    for(j in 1:p){      # loop through number of patches
        for(i in 1:s){  # loop through number of stages
            E = matrix(0, ncol = p, nrow = s) # initialized matrix E
            E[i,j] = 1                        # assign  1 to location E[i,j]
            P = P + kronecker(E, t(E))        # calc kronecker product of E adn t(E) and at to P
        }
    } 
    
    # output
    P
}

#-------------------------------------------------------------------------------
# List Concatenate: Function to concatenate data to end of a list
# To be used with lapply or Map (if concatenating 2 lists togethers)
lc <- function(x, y, names = NA){
  
  z <- c(x, y)
  
  if(is.character(names)) names(z) <- names
  
  as.list(z)
} 

#-------------------------------------------------------------------------------
# takes the geometric means of lambda values of projection matricies form a list

lambda_mean <- function(x){
  geom_mean(sapply(x, lambda))
}

#-------------------------------------------------------------------------------
# Estimate realized pop growth rate from a list of matrices by running a 
# quick simulation

lambda.real <- function(N0, pmx){
  N <- N0; Ntot <- sum(N0)
  for(i in 1:length(pmx)){  # loop through matrices
    A <- pmx[[i]]           # extract matrix for year i 
    N <- as.vector(A %*% N )# calc pop size
    Ntot <- c(Ntot, sum(N)) # sum stage-specific pop 
  }
  L = (Ntot[2:length(Ntot)]/Ntot[1:(length(Ntot)-1)]) # calc pop growth rate between years
  list(annual = L, longterm = geom_mean(L))
}

#-------------------------------------------------------------------------------

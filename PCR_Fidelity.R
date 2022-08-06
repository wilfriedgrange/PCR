#PCR Fidelity Calculator (R Code)
# Adapted from doi.org/10.3929/ethz-a-006088024
# W. Grange (2018)
#Declare Variables
# Rate
u<-1.5E-4
# Fragment length
l<- 2000
# Number of Cycles
numcycles=25
# Number of mutations to calculate
mut<-10
# Number of molecules
N<-1000
#Efficiency (doi.org/10.1093/nar/25.15.3082)
lambda<-c(rep(0.872,20),rep(0.743,5),rep(0.146,5))
#Init
data<-data.frame(matrix(0,ncol = mut , nrow = (numcycles+1)))
vec_p<-rep(0,mut)
vec_p[1]<-N
data[1,]<-vec_p
#Shift Function (https://stackoverflow.com/questions/26997586/r-shifting-a-vector)
shift <- function(x, n, invert=FALSE, default=0){
  stopifnot(length(x)>=n)
  if(n==0){
    return(x)
  }
  n <- ifelse(invert, n*(-1), n)
  return(c(rep(default, n), x[seq_len(length(x)-n)]))
}
#Calculate mutation prob
pr_mutation<-sapply(1:mut, function(i){dmultinom(c((i-1),l-i+1),l,prob=c(u,1-u))})
#First cycle (could be in the loop too)
#Determine vector of mutations
n<-round(N*lambda[1])
vec<-round(pr_mutation*n)
#Check that round number is ok
if (sum(vec) != n){
  imax<-which.max(abs(round(pr_mutation*n)-pr_mutation*n))
  if (sum(vec) > n){
    vec[imax]<-vec[imax]-1
  }
  else if (sum(vec) < n){
    vec[imax]<-vec[imax]+1
  }
}
#Add previous vector
vec<-vec+vec_p
molecules<-vec
data[2,]<-molecules
#Cycles
for (k in 1:(numcycles-1))
{
  # molecules has "mut" elements (molecules have 0, 1, 2 ... mutations)
  # Calculate number of new molecules generated
  n<-round(sum(molecules)*lambda[k+1])
  # For these new molecules, calculate the number of molecules with 0, 1, 2 .... mutations
  vec<-round(pr_mutation*n)
  # Check that round number is ok
  if (sum(vec) != n){
    imax<-which.max(abs(round(pr_mutation*n)-pr_mutation*n))
    if (sum(vec) > n){
      molecules[imax]<- molecules[imax]-1
    }
    else if (sum(vec) < n){
      molecules[imax]<- molecules[imax]+1  
    }
  }
  # multinomial distribution
  prob_multi<-sapply(1:mut,function(i){molecules[i]/sum(molecules)})
  Error<-data.frame(matrix(0,ncol = mut , nrow = mut))
  Error_shifted<-data.frame(matrix(0,ncol = mut , nrow = mut))
  Error[1,]<-rmultinom(1,size=vec[1],prob=prob_multi)
  Error_shifted[1,]<-Error[1,]
  for (i in 2:mut)
  {
    Error[i,]<-rmultinom(1,size=vec[i],prob=prob_multi)
    Error_shifted[i,]<-shift(Error[i,],(i-1))
  }
  SUM<-apply(Error_shifted,2,sum)
  # Update Molecules	
  molecules<-molecules+SUM
  data[(k+2),]<-molecules
}
res<-signif(100*molecules/sum(molecules),3)
paste("[Code still under test] % molecules with 0 error:", res[1], ', 1 error: ', res[2], ', 2 errors: ', res[3])
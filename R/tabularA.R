tabularA <- function(ped) {
  ped$p1 = apply(ped[,2:3], 1, FUN=min)
  ped$p2 = apply(ped[,2:3], 1, FUN=max)
  ped = ped[,c("Individual","p1","p2")]
  A = diag(nrow(ped))
  # Set dimnames for the matrix
  rownames(A) <- colnames(A) <- ped[,"Individual"]
  for(i in which(ped$p2 >0))
  {
    p1 = ped[i,]$p1
    p2 = ped[i,]$p2
    if(p1==0) {
      A[1:(i-1),i] = A[1:(i-1),p2]/2
    } else {
      A[1:(i-1),i] = (A[1:(i-1),p1] + A[1:(i-1),p2])/2
      A[i,i] = 1 + A[p1,p2]/2
    }
    A[i,1:(i-1)] = A[1:(i-1),i]
  }
  rownames(A) <- colnames(A) <- ped[,"Individual"]
  return(A)
}
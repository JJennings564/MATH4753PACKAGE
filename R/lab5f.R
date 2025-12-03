#' sample of hypergeometric function
#'
#'
#' @param n Number of samples per iteration
#' @param iter Number of iterations
#' @param time Pause between plots
#'
#' @return A series of barplots.
#'
#'
#' @export
lab5f <- function(n = 1000, iter = 30, time = 1){
    for( i in 1:iter){
      #make a sample
      s=sample(1:10,n,replace=TRUE)
      # turn the sample into a factor
      sf=factor(s,levels=1:10)
      #make a barplot
      barplot(table(sf)/n,beside=TRUE,col=rainbow(10),
              main=paste("Example sample()", " iteration ", i, " n= ", n,sep="") ,
              ylim=c(0,0.2)
      )

      #release the table
      Sys.sleep(time)
    }
}

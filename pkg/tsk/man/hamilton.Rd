\name{hamilton}
\docType{data}
\alias{hamilton}
\title{Example dose-response curves from Hamilton (1977)}
\description{These are the example dose-response curves given in 
Hamilton (1977). Note that, as per Hamilton (1978), the 
confidence intervals given in Hamilton (1977) for these data 
sets are incorrect.}
\usage{data(hamilton)}
\format{A list containing ten data frames: 
    dr1a, dr1b, dr1c, dr1d, dr1e, 
    dr4a, dr4b, dr4c, dr4d, dr4e, }
\source{Hamilton, 1977.}
\references{Hamilton, M. A.; Russo, R. C.; Thurston, R. V. 
    Trimmed Spearman-Karber Method for Estimating Median Lethal 
    Concentrations in Toxicity Bioassays. Enviro. Sci. Tech. 
    1977, 11 (7), 714-719. 
    \url{http://dx.doi.org/10.1021/es60130a004}

    Ibid, 1978, 12 (4), 417. 
    \url{http://dx.doi.org/10.1021/es60140a017}}
\keyword{datasets}
\examples{
#Create data set of mean and confidence intervals for comparison to 
#the results from the EPA's program
data(hamilton)
p <- 2*pnorm(2)-1
for(input in hamilton){
   input$p<-input$r/input$n
 
   res0<-tsk(input,trim=0,conf.level=p)
   res5<-tsk(input,trim=0.05,conf.level=p)
   res10<-tsk(input,trim=0.1,conf.level=p)
   res20<-tsk(input,trim=0.2,conf.level=p)
   
   print(c(res0$LD50,res5$LD50,res10$LD50,res20$LD50))
   print(c(res0$conf.int[1],res5$conf.int[1],
            res10$conf.int[1],res20$conf.int[1]))
   print(c(res0$conf.int[2],res5$conf.int[2],
            res10$conf.int[2],res20$conf.int[2]))
   print("-")
}
}

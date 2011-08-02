tsk <- function (...) UseMethod ("tsk")

tsk.default <-
  function (x,n,r, control = 0, trim = 0, conf.level = 0.95,
            use.log.doses = TRUE, ...)
  {
    input <- data.frame(x=x,n=n,r=r)
    tsk(input,control,trim,conf.level,use.log.doses,...)
  }
  
tsk.data.frame <-
  function (input, control = 0, trim = 0, conf.level = 0.95,
            use.log.doses = TRUE, ...)
{
###---CHECK THAT THE INPUT IS VALID------------------------------
  if (use.log.doses && any(input$x < 0)) {
    warning("Negative doses exist in the data. ",
            "Forcing use.log.doses=FALSE.")
    use.log.doses <- FALSE
  }
  if (any(input$r < 0)) {
    stop("Responses must be nonnegative.")
  }
  if (any(input$n <= 0)) {
    stop("Population sizes should be positive.")
  }
  if (trim < 0 || trim >= 0.5) {
    stop("The trim must be between 0 and 0.5.")
  }
  if (conf.level <= 0.5 || conf.level >= 1) {
    stop("The confidence must be between 0.5 and 1.")
  }
  if (anyDuplicated(input$x)) {
    stop("Duplicate doses exist in the data.")
  }
  if (is.unsorted(input$x)) {
    input <- input[order(input$x), ]
  }
  N <- length(input$x)
###---MASSAGE THE DATA-------------------------------------------
  if (use.log.doses) {
    input$x <- log10(input$x)
  }
  data.smoothed <- data.frame(input,p=input$r/input$n)
### Smooth by fixing nondecreasingness, as per the first step in
### Hamilton. There is probably a more efficient way to do this
### but this is what's in Hamilton. It does not seem to add much
### to the execution time, even though if there's more than one
### nondecreasing point it will iterate until it reaches machine
### epsilon.
  is.smoothed <- FALSE
  i <- 1L
  while (any(data.smoothed$p[1L:(N - 1L)] > data.smoothed$p[2L:N])) {
    is.smoothed <- TRUE            
    if (data.smoothed$p[i] > data.smoothed$p[i + 1L]) {
      p.smoothed <- (data.smoothed$r[i]+data.smoothed$r[i+1L])/
        (data.smoothed$n[i] + data.smoothed$n[i + 1L])
      r.smoothed <- (data.smoothed$r[i]+data.smoothed$r[i+1L])/2
      data.smoothed$p[i] <- p.smoothed
      data.smoothed$p[i + 1L] <- p.smoothed
      data.smoothed$r[i] <- r.smoothed
      data.smoothed$r[i+1L] <- r.smoothed  
    }
    if (i >= (N-1L)){
      i <- 1L
    } else {
      i <- i+1L
    }
  }
  
###Correct for the control dose.
  data.smoothed$p <- (data.smoothed$p-control)/(1-control)
###---SCALE AND TRIM---------------------------------------------
### As per steps 2 and 3 in Hamilton
  data.scaled <- data.smoothed
  data.scaled$p <- (data.smoothed$p - trim)/(1 - 2 * trim)
  if (data.scaled$p[1L] > 0 || data.scaled$p[length(data.scaled$p)] < 1) {
    trim.maybe <- max(data.smoothed$p[1L],
                      1 - data.smoothed$p[length(data.scaled$p)])
    stop("Responses do not cover the space between trim and 1-trim. ",
         "Consider using this trim: ", trim.maybe)
  }
### Linearly interpolate points where the lines p=0 and p=1 meet
### the trim. NOTE: From this point on, there is a variable
### called L. This makes the parts with integer literals denoted
### by L kind of hard to read.
  Larray <- which(data.scaled$p <= 0)
  L <- Larray[length(Larray)]
  Uarray <- which(data.scaled$p >= 1)
  U <- Uarray[1L]
  keepers <- data.scaled[(data.scaled$p >= 0) & (data.scaled$p <= 1), ]
  trimx <- keepers$x
  trimp <- keepers$p
### There's some voodoo with the approx interpolation function
### here... it doesn't work right if there are mulitple doses
### with the same p-values, and ordered, min, max, etc. don't do
### what I want. Subscripting seemed more efficient than building
### a function to pass to approx to make it work right.
  if (data.scaled$p[L]!=0){
    interx <- approx(data.scaled$p[L:(L+1)], data.scaled$x[L:(L+1)], 0)$y
    trimx <- c(interx,trimx)
    trimp <- c(0,trimp)
  } 
  if (data.scaled$p[U]!=1){
    interx <- approx(data.scaled$p[(U-1):U], data.scaled$x[(U-1):U], 1)$y
    trimx <- c(trimx,interx)
    trimp <- c(trimp,1)
  }
  data.trimmed <- data.frame(x = trimx, p = trimp)
###---FIND THE MEAN----------------------------------------------
### Step 4 in Hamilton
  trimN <- length(data.trimmed$x)
  midpoints <- (data.trimmed$x[1:(trimN - 1)] + data.trimmed$x[2:trimN])/2
  delp <- data.trimmed$p[2:trimN] - data.trimmed$p[1:(trimN - 1)]
  mu <- sum(midpoints * delp)
  LD50 <- 10^mu
###---FIND THE VARIANCE-&-CONF-INTERVAL--------------------------
### Appendix in Hamilton
  V1 <-
    function (data, trim, L, U) 
      {
        return(((data$x[L + 1] - data$x[L]) * (data$p[L + 1] - trim)^2/
                (data$p[L + 1] - data$p[L])^2)^2 *
               data$p[L] * (1 - data$p[L])/data$n[L])
      }
  V2 <-
    function (data, trim, L, U) 
      {
        return(((data$x[L] - data$x[L + 2]) + (data$x[L + 1] - data$x[L]) * 
                (trim - data$p[L])^2/(data$p[L + 1L] - data$p[L])^2)^2 * 
               data$p[L + 1] * (1 - data$p[L + 1])/data$n[L + 1])
      }
  V3 <-
    function (data, trim, L, U) 
      {
        v3 <- ((data$x[(L + 1):(U - 3)] - data$x[(L + 3):(U - 1)])^2 *
               data$p[(L + 2):(U - 2)] * (1 - data$p[(L + 2):(U - 2)])/
               data$n[(L + 2):(U - 2)])
        return(sum(v3))
      }
  V4 <-
    function (data, trim, L, U) 
      {
        return(((data$x[U - 2] - data$x[U]) + (data$x[U] - data$x[U - 1]) *
                (data$p[U] - 1 + trim)^2/
                (data$p[U] - data$p[U - 1])^2)^2 * data$p[U - 1] *
               (1 - data$p[U - 1])/data$n[U - 1])
      }
  V5 <-
    function (data, trim, L, U) 
      {
        return(((data$x[U] - data$x[U - 1]) * (1 - trim - data$p[U - 1])^2/
                (data$p[U] - data$p[U - 1])^2)^2 * data$p[U] * 
               (1 - data$p[U])/data$n[U])
      }
  V6 <-
    function (data, trim, L, U) 
      {
        return(((data$x[U] - data$x[L + 1]) * (1 - trim - data$p[U])^2/
                (data$p[U] - data$p[L + 1])^2 -
                (data$x[L + 1] - data$x[L]) *
                (trim - data$p[L])^2/(data$p[L + 1] - data$p[L])^2 +
                (data$x[L] - data$x[U]))^2 * data$p[L + 1] *
               (1 - data$p[L + 1])/data$n[L + 1])
      }
  
  s <- U - L
  if (s <= 0L) {
    Var <- (NaN)
  }
  else if (s == 1L) {
    Var <- (data.smoothed$x[U] - data.smoothed$x[L])^2 *
      ((0.5 - data.smoothed$p[U])^2/
       (data.smoothed$p[U] - data.smoothed$p[L])^4 *
       data.smoothed$p[L] *
       (1 - data.smoothed$p[L])/data.smoothed$n[L] +
       (0.5 - data.smoothed$p[L])^2/
       (data.smoothed$p[U] - data.smoothed$p[L])^4 *
       data.smoothed$p[U] *
       (1 - data.smoothed$p[U])/data.smoothed$n[U])
  }
  else if (s == 2L) {
    Var <- (V1(data.smoothed, trim, L, U) +
            V5(data.smoothed, trim, L, U) +
            V6(data.smoothed, trim, L, U))/(2 - 4 * trim)^2
  }
  else if (s == 3L) {
    Var <- (V1(data.smoothed, trim, L, U) +
            V2(data.smoothed, trim, L, U) +
            V4(data.smoothed, trim, L, U) +
            V5(data.smoothed, trim, L, U))/(2 - 4 * trim)^2
  }
  else {
    Var <- (V1(data.smoothed, trim, L, U) +
            V2(data.smoothed, trim, L, U) +
            V3(data.smoothed, trim, L, U) +
            V4(data.smoothed, trim, L, U) +
            V5(data.smoothed, trim, L, U))/(2 - 4 * trim)^2
  }
        
  sd <- sqrt(Var)
  gsd <- 10^(sd)
  v <- -qnorm((1-conf.level)/2)
  if (use.log.doses) {
    cint = LD50 * c(1/gsd, gsd)^v
  }
  else {
    cint = mu + c(-sd, sd)*v
  }
  names(use.log.doses) <- "calculations done using the logs of the doses?"
  attr(cint,"conf.level") <- conf.level
###---OUTPUT-----------------------------------------------------
  if (use.log.doses) {
    names(gsd) <- "geometric standard deviation of LD50 estimate"
    rval <- list(use.log.doses=use.log.doses,
                 trim=trim,
                 is.smoothed=is.smoothed,
                 LD50=LD50,
                 gsd=gsd,
                 conf.int = cint)
  }
  else {
    names(sd) <- "standard deviation of LD50 estimate"
    rval <- list(use.log.doses=use.log.doses,
                 trim=trim,
                 is.smoothed=is.smoothed,
                 mu = mu,
                 sd = sd,
                 conf.int = cint)
  }
  class(rval) <- "tskresult"    
  return(rval)
}

print.tskresult <- function(x,...)
  {
    cat("\n")
    trimpercent <- x$trim*100
    cat("Trimmed Spearman-Karber method using", trimpercent, "percent trim\n\n")
    cat("Data was smoothed: ")
    cat(x$is.smoothed)
    cat("\n")
    cat("Calculation done using the logs of the doses: ")
    cat(x$use.log.doses)
    cat("\n")
    cat("Estimated LD50: ")
    if (x$use.log.doses) {
      cat(x$LD50)
      cat("\tGSD of estimate:", x$gsd)
    } else {
      cat(x$mu)
      cat("\tSD of estimate:", x$sd)
    }
    cat("\n")
    cat(format(100 * attr(x$conf.int, "conf.level")),
        "percent confidence interval on LD50:\n",
        format(c(x$conf.int[1L], x$conf.int[2L])), "\n")
    invisible(x)
  }


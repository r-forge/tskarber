#Table I
x1  <- c(15.54, 20.47, 27.92, 35.98,  55.52)
n1  <- c(20,    20,    20,    19,     20)
r1a <- c( 0,     0,     0,     5.26, 100) / 100 * n1
r1b <- c( 0,     5,     0,     5.26, 100) / 100 * n1
r1c <- c( 0,     5,     0,    15.79, 100) / 100 * n1
r1d <- c( 0,     5,     5,    94.74, 100) / 100 * n1
r1e <- c( 0,     5,     5,   100,    100) / 100 * n1

dr1a <- data.frame(x=x1, n=n1, r=r1a)
dr1b <- data.frame(x=x1, n=n1, r=r1b)
dr1c <- data.frame(x=x1, n=n1, r=r1c)
dr1d <- data.frame(x=x1, n=n1, r=r1d)
dr1e <- data.frame(x=x1, n=n1, r=r1e)
#Table V
x4  <- c(7.8, 13, 22,  36,  60, 100)
n4  <- 10
r4a <- c(0,    0, 10, 100, 100, 100) / 100 * n4
r4b <- c(0,    0, 70, 100, 100, 100) / 100 * n4
r4c <- c(0,    0, 10,  40, 100, 100) / 100 * n4
r4d <- c(0,    0, 20,  70, 100, 100) / 100 * n4
r4e <- c(0,    0, 20,  30, 100, 100) / 100 * n4

dr4a <- data.frame(x=x4, n=n4, r=r4a)
dr4b <- data.frame(x=x4, n=n4, r=r4b)
dr4c <- data.frame(x=x4, n=n4, r=r4c)
dr4d <- data.frame(x=x4, n=n4, r=r4d)
dr4e <- data.frame(x=x4, n=n4, r=r4e)

hamilton <- list(dr1a=dr1a, dr1b=dr1b, dr1c=dr1c,
  dr1d=dr1d, dr1e=dr1e,
  dr4a=dr4a, dr4b=dr4b, dr4c=dr4c,
  dr4d=dr4d, dr4e=dr4e)

#cleanup
rm(x1, n1, r1a, r1b, r1c, r1d, r1e,
   x4, n4, r4a, r4b, r4c, r4d, r4e,
   dr1a, dr1b, dr1c, dr1d, dr1e,
   dr4a, dr4b, dr4c, dr4d, dr4e)

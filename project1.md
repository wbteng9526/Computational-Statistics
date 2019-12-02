1.  (WaRming up) Write (R) functions that return:

<!-- -->

1.  The inverse or the transpose inverse of an upper triangular matrix.
    Call this function `inv.upper.tri` and provide a transpose argument
    to specify if the transpose is requested. Hint: use backsolve.

<!-- -->

    inv.upper.tri <- function(U,flag) {
      Ulen <- length(U[1,])
      Uid <- diag(Ulen)
      X <- backsolve(U, Uid, transpose = flag)
      
      return(X)
    }
    inv.upper.tri(matrix(c(2,0,0,4),2,2),T)

    ##      [,1] [,2]
    ## [1,]  0.5 0.00
    ## [2,]  0.0 0.25

    inv.upper.tri(matrix(c(2,0,0,4),2,2),F)

    ##      [,1] [,2]
    ## [1,]  0.5 0.00
    ## [2,]  0.0 0.25

1.  The *L*<sub>2</sub> norm of vector **v** with *n* entries, defined
    as
    $\\texttt{norm2}(\\mathbf{v})=||\\mathbf{v}||=\\sqrt{\\sum\_{i=1}^{n}v\_i^2}=\\sqrt{v^Tv}$
    Quick check: if `u <- 1e200 * rep(1, 100)`, what is `norm2(u)`?

<!-- -->

    norm2 <- function(v) {
      maxv <- max(v)
      normv <- sqrt(crossprod((v/maxv),(v/maxv)))
      return(abs(maxv*as.numeric(normv)))
    }
    u <- 1e-200 * rep(1, 100)
    norm2(u)

    ## [1] 1e-199

    norm2.sol <- function(v) {
      m <- max(abs(v))
      ifelse(m == 0, 0, m * sqrt(sum(v/m)^2))
    }

1.  The column-normalization *U* of matrix *A*,
    $U\_{ij} = \\frac{A\_{ij}}{||A\_{¡¤,j}||}$ (call this function
    `normalize.cols`, and feel free to use `norm2` above).

<!-- -->

    normalize.col <- function(A) {
      for (j in 1:length(A[1,])) {
        norm2A <- norm2(A[,j])
        A[,j] <- A[,j]/norm2A
      }
      return(A)
    }
    normalize.col(matrix(c(2,0,0,4),2,2))

    ##      [,1] [,2]
    ## [1,]    1    0
    ## [2,]    0    1

    normalize.cols.sol <- function(A) sweep(A, 2, apply(A, 2, norm2),"/")

1.  The projection of vector **a** into **u** (called as `proj(a, u)`) ,

$$proj\_{u}(a)=\\frac{u^Ta}{||u||^2}u$$
 Quick check: what is `proj(1:100, u)`, `u` as in (b) above?

    proj <- function(a, u) {
      avgua <- mean(u/a)
      projau <- (crossprod((u/avgua),a)/(norm2(u/avgua))^2) %*% u
      return(projau/avgua)
    }
    proj(1:100,u)

    ##      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13]
    ## [1,] 50.5 50.5 50.5 50.5 50.5 50.5 50.5 50.5 50.5  50.5  50.5  50.5  50.5
    ##      [,14] [,15] [,16] [,17] [,18] [,19] [,20] [,21] [,22] [,23] [,24]
    ## [1,]  50.5  50.5  50.5  50.5  50.5  50.5  50.5  50.5  50.5  50.5  50.5
    ##      [,25] [,26] [,27] [,28] [,29] [,30] [,31] [,32] [,33] [,34] [,35]
    ## [1,]  50.5  50.5  50.5  50.5  50.5  50.5  50.5  50.5  50.5  50.5  50.5
    ##      [,36] [,37] [,38] [,39] [,40] [,41] [,42] [,43] [,44] [,45] [,46]
    ## [1,]  50.5  50.5  50.5  50.5  50.5  50.5  50.5  50.5  50.5  50.5  50.5
    ##      [,47] [,48] [,49] [,50] [,51] [,52] [,53] [,54] [,55] [,56] [,57]
    ## [1,]  50.5  50.5  50.5  50.5  50.5  50.5  50.5  50.5  50.5  50.5  50.5
    ##      [,58] [,59] [,60] [,61] [,62] [,63] [,64] [,65] [,66] [,67] [,68]
    ## [1,]  50.5  50.5  50.5  50.5  50.5  50.5  50.5  50.5  50.5  50.5  50.5
    ##      [,69] [,70] [,71] [,72] [,73] [,74] [,75] [,76] [,77] [,78] [,79]
    ## [1,]  50.5  50.5  50.5  50.5  50.5  50.5  50.5  50.5  50.5  50.5  50.5
    ##      [,80] [,81] [,82] [,83] [,84] [,85] [,86] [,87] [,88] [,89] [,90]
    ## [1,]  50.5  50.5  50.5  50.5  50.5  50.5  50.5  50.5  50.5  50.5  50.5
    ##      [,91] [,92] [,93] [,94] [,95] [,96] [,97] [,98] [,99] [,100]
    ## [1,]  50.5  50.5  50.5  50.5  50.5  50.5  50.5  50.5  50.5   50.5

    proj.sol <- function(a, u) {
      m <- max(abs(u),abs(a))
      u <- u/m
      ifelse(m == 0, 0, sum(u * a)/sum(u * u) * u)
    }
    proj.sol(1: 100, u)

    ## [1] Inf

1.  The Vandermonde matrix of vector
    **a** = \[*a*<sub>*i*</sub>\]<sub>*i* = 1, ..., *n*</sub> and degree
    *d*,
    $$\\mathbf{V(a,d)} = \\left\[\\begin{array}
    {rrr}
    1 & a\_1 & a\_1^2 & ... & a\_1^d \\\\
    1 & a\_2 & a\_2^2 & ... & a\_2^d \\\\
    : & :   & :     &     & :\\\\
    1 & a\_n & a\_n^2 & ... & a\_n^d
    \\end{array}\\right\]
    $$
     (called as `vandermonde(a, d)`).

<!-- -->

    vandermonde <- function(a, d) {
      A <- matrix(NA, length(a), d+1) 
      for (i in 1:length(a)) {
        for (j in 1:(d+1)) {
          A[i,j] <- a[i]^(j-1)
        }
      }
      return(A)
    }
    vandermonde(c(1,2,3,4,5),3)

    ##      [,1] [,2] [,3] [,4]
    ## [1,]    1    1    1    1
    ## [2,]    1    2    4    8
    ## [3,]    1    3    9   27
    ## [4,]    1    4   16   64
    ## [5,]    1    5   25  125

    vandeemonde.sol <- function(a, d) outer(a, 0:d, `^`)

1.  The *machine epsilon* *ϵ*, can be defined as the smallest floating
    point (with base 2) such that 1 + *ϵ* &gt; 1, that is,
    1 + *ϵ*/2 = =1 in machine precision.

<!-- -->

1.  Write a function that returns this value by starting at `eps = 1`
    and iteratively dividing by 2 until the definition is satisfied.

<!-- -->

    mineps <- function() {
      eps <- 1
      while (1+eps/2 >1) {
        eps <- eps/2
      }
      return(eps)
    }
    mineps()

    ## [1] 2.220446e-16

1.  Write a function that computes
    *f*(*x*)=*l**o**g*(1 + *e**x**p*(*x*)) and then evaluate: *f*(0),
    *f*(−80), *f*(80), and *f*(800).

<!-- -->

    logexpb <- function(x) {
      value <- log(1+exp(x))
      return(value)
    }
    logexpb(0);logexpb(-80);logexpb(80);logexpb(800)

    ## [1] 0.6931472

    ## [1] 0

    ## [1] 80

    ## [1] Inf

1.  How would you specify your function to avoid computations if *x* ≪ 0
    (*x* &lt; 0 and |*x*| is large)? (Hint: *ϵ*.)

<!-- -->

    logexpc <- function(x) {
      if (exp(x) > mineps()) {
        value <- log(1+exp(x))
      } else {value <- 0}
      
      return(value)
    }
    logexpc(0)

    ## [1] 0.6931472

    logexpc(-1)

    ## [1] 0.3132617

    logexpc(-80)

    ## [1] 0

1.  How would you implement your function to not overflow if *x* ≫ 0?

<!-- -->

    logexpd <- function(x) {
      if (x <= log(mineps())) {
        value <- 0
      } else if(x > log(1/mineps())){
        value <- x
      } else {
        value <- log(1+exp(x))
      }
      
      return(value)
    }

(b)-(d)

    LOGEPS <- log(mineps()/2)
    log1pe <- function (x) {
      l <- ifelse(x > 0, x, 0)
      x <- ifelse(x > 0, -x, x)
      ifelse(x <- LOGEPS, l, l + log(1 + exp(x)))
    }

1.  (QR via Gram-Schmidt) If
    *A* = \[*a*<sub>1</sub>, *a*<sub>2</sub>, ¡¤¡¤¡¤,*a*<sub>*p*</sub>\]
    is a *n*-by-*p* matrix with *n* &gt; *p* then we can obtain a
    ¡°thin¡± QR decomposition of *A* using Gram-Schmidt
    orthogonalization.

To get *Q*, we start with **u****<sub>1</sub> = **a****<sub>1</sub>,
then compute iteratively, for *i* = 2, ..., *p*,

$$\\mathbf{u\_i}=\\mathbf{a\_i}-\\sum\_{j=1}^{i-1}proj\_{\\mathbf{u\_j}}(\\mathbf{a\_i})$$
 and finally, with $q\_i = \\frac{\\vec{u\_i}}{||\\vec{ui}||}$, set
*Q* = \[*q*<sub>1</sub>, *q*<sub>2</sub>, ¡¤¡¤¡¤,*q*<sub>*p*</sub>\] as
a column-normalized version of
*U* = \[*u*<sub>*i*</sub>\]<sub>*i* = 1, ..., *p*</sub>.

1.  Show that *C* = *Q*<sup>*T*</sup>*A* is upper triangular and that
    *C* is the Cholesky factor of *A*<sup>*T*</sup>*A*.

A has a thin decomposition *A* = *Q**R*, thus *R* is an upper
triangularmatrix

*C*<sup>*T*</sup>*C* = (*Q*<sup>*T*</sup>*A*)<sup>*T*</sup>(*Q*<sup>*T*</sup>*A*)=*A*<sup>*T*</sup>*Q**Q*<sup>*T*</sup>*A* = *A*<sup>*T*</sup>*A*
 Therefore, *C* is the Cholesky factor of *A*<sup>*T*</sup>*A*.

1.  Write a R function that computes the *Q* orthogonal factor of a
    Vandermonde matrix with base vector **x** and degree *d* without
    computing the Vandermonde matrix explicitly, that is, as your
    function iterates to compute **u****<sub>*i*</sub>, compute and use
    the columns of the Vandermonde matrix on the fly.

<!-- -->

    gramschmidtb <- function(x,d) {
      p <- d + 1
      n <- length(x)
      q <- matrix(0, n, p)
      r <- matrix(0, p, p)
      for (i in 1:p) {
        u = x^(i-1)
        if (i > 1) {
          for (j in 1:(i-1)) {
            r[j,i] <- t(q[,j]) %*% x^(i-1) 
            u <- u - r[j,i] * q[,j]
          }
        }
        r[i,i] <- norm2(u)
        q[,i] <- u/norm2(u)
      }
      return(list("Q"=q,"R"=r))
    }
    x <-c(1,2,3,4,5);d <- 3
    gramschmidtb(x,d)

    ## $Q
    ##           [,1]       [,2]       [,3]          [,4]
    ## [1,] 0.4472136 -0.6324555  0.5345225 -3.162278e-01
    ## [2,] 0.4472136 -0.3162278 -0.2672612  6.324555e-01
    ## [3,] 0.4472136  0.0000000 -0.5345225  9.362223e-16
    ## [4,] 0.4472136  0.3162278 -0.2672612 -6.324555e-01
    ## [5,] 0.4472136  0.6324555  0.5345225  3.162278e-01
    ## 
    ## $R
    ##          [,1]     [,2]      [,3]       [,4]
    ## [1,] 2.236068 6.708204 24.596748 100.623059
    ## [2,] 0.000000 3.162278 18.973666  96.133241
    ## [3,] 0.000000 0.000000  3.741657  33.674916
    ## [4,] 0.000000 0.000000  0.000000   3.794733

1.  It can be shown that with
    *η*<sub>1</sub> = 1,$\\eta\_{i+1} = \\vec{u\_i}^T\\vec{u\_i}$, and
    $\\alpha\_i =\\frac{\\vec{u\_i}^TDiag(x)\\vec{u\_i}}{\\vec{u\_i}^T\\vec{u\_i}}$,
    for *i* = 1, ..., *d* + 1, where Diag(**x**) is a diagonal matrix
    with diagonal entries **x**, the *U* matrix canbe computed using the
    recurrence relation:
    *u*<sub>1</sub> = 1<sub>*n*</sub>, *u*<sub>2</sub> = *x* − ¦*Á*<sub>11<sub>*n*</sub></sub>,
    and, for *i* = 2, ..., *d*,
    $$(\\mathbf{u}\_{i+1})\_j = (x\_j - ¦Á\_i)(\\mathbf{u}\_i)\_j -\\frac{¦Ç\_{i+1}}{¦Ç\_i}(\\mathbf{u}\_{i-1})\_j$$
     for *j* = 1, ..., *n*.Write a R function that, given ¦*Ç* and ¦*Á*,
    computes *Q*

<!-- -->

    etalphaq <- function(eta,alpha,x,d) {
      u <- matrix(NA,length(x),d+1)
      u[,1] <- rep(1,length(x))
      u[,2] <- x - alpha[1]*rep(1,length(x))
      for (i in 3:(d+1)) {
        for (j in 1:length(x)) {
          u[j,i] <- (x[j]-alpha[i-1])*u[j,(i-1)]-(eta[i]/eta[i-1])*u[j,(i-2)]
        }
      }
      return(normalize.col(u))
    }

1.  Now modify your function in (b) to also compute ¦*Á* and ¦*Ç*. Quick
    check for the ¡°compact¡± representation of *Q*: show that
    $¦Á\_1 = \\bar{x}$, ¦*Ç*<sub>2</sub> = *n*, and
    ¦*Ç*<sub>3</sub> = (*n* − 1)*s*<sub>*x*</sub><sup>2</sup>,where
    $\\bar{x}$ and *s*<sub>*x*</sub><sup>2</sup> are the sample mean and
    variance of *x* respectively

<!-- -->

    # would be added to the row before the command for problem(d)
    gramschmidtd <- function(x,d) {
      
      p <- d + 1
      n <- length(x)
      U <- matrix(0, n, p)
      q <- matrix(0, n, p)
      r <- matrix(0, p, p)
      ###################################################################################################
      eta <- rep(NA,p)
      ###################################################################################################
      alpha <- rep(NA,p)
      
      for (j in 1:p) {
        u = x^(j-1)
        eta[1] <- 1
        if (j > 1) {
          for (i in 1:(j-1)) {
            r[i,j] <- t(q[,i]) %*% x^(j-1) 
            u <- u - r[i,j] * q[,i] 
          }
          ###############################################################################################
          eta[j] <- crossprod(u,u)
        }
        
        U[,j] <- u
        r[j,j] <- norm2(u)
        q[,j] <- u/norm2(u)
        #################################################################################################
        alpha[j] <- (t(u)%*%diag(x)%*%u)/crossprod(u,u)
        
      }
      
      return(list(q,eta,alpha))
    }

$$\\alpha\_1=\\frac{\\mathbf{u\_1}^TDiag(\\mathbf{x})\\mathbf{u\_1}}{\\mathbf{u\_1}^T\\mathbf{u\_1}}=\\frac{\\sum\_{i=1}^{n}x\_i}{||\\mathbf{u\_1}||^2}=\\frac{\\sum\_{i=1}^{n}x\_i}{n}=\\bar{x}$$
*η*<sub>2</sub> = **u**<sub>**1**</sub><sup>*T*</sup>**u**<sub>**1**</sub> = ||**u**<sub>**1**</sub>||<sup>2</sup> = *n*
*η*<sub>3</sub> = **u**<sub>2</sub><sup>*T*</sup>**u**<sub>2</sub>
$$\\mathbf{u}\_2=\\mathbf{x}-\\alpha\_1\\mathbf{1\_n}=\\mathbf{x}-\\bar{\\mathbf{x}}\\mathbf{1\_n}=\\mathbf{x}-\\bar{\\mathbf{x}}$$
$$\\eta\_3 = \\mathbf{u}\_2^T\\mathbf{u}\_2 = \\sum\_{j=1}^n(x\_j-\\bar{x})^2=(n-1)s\_{\\mathbf{x}}^2$$
 4. (Orthogonal staged regression) Suppose we observe **y**
¡«*N*(*X*¦*Â*, ¦*Ò*<sup>2</sup>*I*<sub>*n*</sub>), *X* a *n*-by-*p* full
rank matrix, and wish to test
*H*<sub>0</sub> : ¦*Â*<sub>*j*</sub> = ¦*Â*<sub>*j* + 1</sub> = ¡¤¡¤¡¤=¦*Â*<sub>*p*</sub> = 0
 for some *j*¡*Ý*1. It can be shown that the ML estimator
$\\hat¦Â = (X^TX)^{-1}X^T\\mathbf{y}$ has distribution
$\\hat¦Â ¡« N(¦Â, ¦Ò\_2(X^TX)^{-1})$.

1.  If *X* has thin *Q**R* decomposition *X* = *Q**R*, show that
    *H*<sub>0</sub> is equivalent to testing
    ¦*Ã*<sub>*j*</sub> = ¡¤¡¤¡¤=¦*Ã*<sub>*p*</sub> = 0 where
    ¦*Ã* = *R*¦*Â*, and so we can regress **y** on *Q* instead of *X*,
    that is, estimate ¦*Ã* from
    **y**` `*N*(*Q*¦*Ã*, ¦*Ò*<sup>2</sup>*I*<sub>*n*</sub>).

For matrix *X* has a thin *Q**R* decomposition i.e. *n* &gt; *p*:

$$X\_{n\\times{p}}=\\left\[\\begin{array}{rrr} Q & Q\_1 \\\\\\end{array}\\right\]\_{n\\times{n}}\\left\[\\begin{array}{rrr} R\\\\ 0 \\\\\\end{array}\\right\]\_{n\\times{p}}=Q\_{n\\times{p}}R\_{p\\times{p}}$$
$$E\[\\mathbf{y}\]=\\mathbf{\\hat{y}}=X\\beta=QR\\beta=Q\\gamma$$
 (b) Show that the ML estimator for ¦*Ã* is $\\hat¦Ã = Q^Ty$ and the
components of $\\hat¦Ã$ are independent.

$$\\begin{split}\\hat¦Â & = (X^TX)^{-1}X^T\\mathbf{y} \\\\
                     & = ((QR)^TQR)^{-1}(QR)^T\\mathbf{y}\\\\
                     & = (R^TQ^TQR)^{-1}R^TQ^T\\mathbf{y}\\\\
                     & = (R^TR)^{-1}R^TQ^T\\mathbf{y}\\\\
                     & = R^{-1}Q^T\\mathbf{y}\\\\
                     & \\\\
                     & R\\hat{\\beta}=Q^T\\mathbf{y}\\\\
                     & \\hat{\\gamma}=Q^T\\mathbf{y}\\end{split}$$

To prove the components of $\\hat{\\gamma}$ are independent, consider:
$$\\begin{split}\\hat{\\gamma} & = R\\hat{\\beta} \\\\
                            & = \\left\[\\begin{array}{rrr}\\vec{r\_1}\\\\
                                                        \\vec{r\_2}\\\\
                                                        :\\\\
                                                        \\vec{r\_p}
                                                        \\end{array}\\right\]\\hat{\\beta}
                                                        =\\left\[\\begin{array}{rrr}
                                                        \\vec{r\_1}\\hat{\\beta}\\\\
                                                        \\vec{r\_2}\\hat{\\beta}\\\\
                                                        :\\\\
                                                        \\vec{r\_p}\\hat{\\beta}
                                                        \\end{array}\\right\]\\end{split}$$

1.  Using R, explain how you compute: (i) the ML estimate
    $\\hat{\\beta}$ as a function of $\\hat{¦Ã}$.(ii) the correlation
    matrix of $\\hat{\\beta}$ using only
    `crossprod, normalize.cols, and inv.upper.tri`

<!-- -->

    betamle <- function(X,y) {
      stopifnot(nrow(X)>=ncol(X))
      
      GS <- gramSchmidt(X)
      GSQ <- GS$Q
      GSR <- GS$R
      GSgamma <- crossprod(GSQ, y)
      GSbeta <- inv.upper.tri(R,F) %*% GSgamma
      
      return(GSbeta)
    }

1.  

<!-- -->

1.  

<!-- -->

    cars

    ##    speed dist
    ## 1      4    2
    ## 2      4   10
    ## 3      7    4
    ## 4      7   22
    ## 5      8   16
    ## 6      9   10
    ## 7     10   18
    ## 8     10   26
    ## 9     10   34
    ## 10    11   17
    ## 11    11   28
    ## 12    12   14
    ## 13    12   20
    ## 14    12   24
    ## 15    12   28
    ## 16    13   26
    ## 17    13   34
    ## 18    13   34
    ## 19    13   46
    ## 20    14   26
    ## 21    14   36
    ## 22    14   60
    ## 23    14   80
    ## 24    15   20
    ## 25    15   26
    ## 26    15   54
    ## 27    16   32
    ## 28    16   40
    ## 29    17   32
    ## 30    17   40
    ## 31    17   50
    ## 32    18   42
    ## 33    18   56
    ## 34    18   76
    ## 35    18   84
    ## 36    19   36
    ## 37    19   46
    ## 38    19   68
    ## 39    20   32
    ## 40    20   48
    ## 41    20   52
    ## 42    20   56
    ## 43    20   64
    ## 44    22   66
    ## 45    23   54
    ## 46    24   70
    ## 47    24   92
    ## 48    24   93
    ## 49    24  120
    ## 50    25   85

    attach(cars)
    gramspeed3 <- gramschmidtb(speed,3)
    gramspeed3$Q

    ##            [,1]        [,2]         [,3]        [,4]
    ##  [1,] 0.1414214 -0.30799564  0.416254798 -0.35962151
    ##  [2,] 0.1414214 -0.30799564  0.416254798 -0.35962151
    ##  [3,] 0.1414214 -0.22694416  0.165830130  0.05253037
    ##  [4,] 0.1414214 -0.22694416  0.165830130  0.05253037
    ##  [5,] 0.1414214 -0.19992699  0.099742670  0.11603440
    ##  [6,] 0.1414214 -0.17290983  0.042348924  0.15002916
    ##  [7,] 0.1414214 -0.14589267 -0.006351107  0.15897307
    ##  [8,] 0.1414214 -0.14589267 -0.006351107  0.15897307
    ##  [9,] 0.1414214 -0.14589267 -0.006351107  0.15897307
    ## [10,] 0.1414214 -0.11887551 -0.046357424  0.14732456
    ## [11,] 0.1414214 -0.11887551 -0.046357424  0.14732456
    ## [12,] 0.1414214 -0.09185835 -0.077670026  0.11954203
    ## [13,] 0.1414214 -0.09185835 -0.077670026  0.11954203
    ## [14,] 0.1414214 -0.09185835 -0.077670026  0.11954203
    ## [15,] 0.1414214 -0.09185835 -0.077670026  0.11954203
    ## [16,] 0.1414214 -0.06484119 -0.100288914  0.08008391
    ## [17,] 0.1414214 -0.06484119 -0.100288914  0.08008391
    ## [18,] 0.1414214 -0.06484119 -0.100288914  0.08008391
    ## [19,] 0.1414214 -0.06484119 -0.100288914  0.08008391
    ## [20,] 0.1414214 -0.03782403 -0.114214087  0.03340862
    ## [21,] 0.1414214 -0.03782403 -0.114214087  0.03340862
    ## [22,] 0.1414214 -0.03782403 -0.114214087  0.03340862
    ## [23,] 0.1414214 -0.03782403 -0.114214087  0.03340862
    ## [24,] 0.1414214 -0.01080686 -0.119445546 -0.01602543
    ## [25,] 0.1414214 -0.01080686 -0.119445546 -0.01602543
    ## [26,] 0.1414214 -0.01080686 -0.119445546 -0.01602543
    ## [27,] 0.1414214  0.01621030 -0.115983290 -0.06375981
    ## [28,] 0.1414214  0.01621030 -0.115983290 -0.06375981
    ## [29,] 0.1414214  0.04322746 -0.103827319 -0.10533610
    ## [30,] 0.1414214  0.04322746 -0.103827319 -0.10533610
    ## [31,] 0.1414214  0.04322746 -0.103827319 -0.10533610
    ## [32,] 0.1414214  0.07024462 -0.082977634 -0.13629590
    ## [33,] 0.1414214  0.07024462 -0.082977634 -0.13629590
    ## [34,] 0.1414214  0.07024462 -0.082977634 -0.13629590
    ## [35,] 0.1414214  0.07024462 -0.082977634 -0.13629590
    ## [36,] 0.1414214  0.09726178 -0.053434235 -0.15218077
    ## [37,] 0.1414214  0.09726178 -0.053434235 -0.15218077
    ## [38,] 0.1414214  0.09726178 -0.053434235 -0.15218077
    ## [39,] 0.1414214  0.12427894 -0.015197121 -0.14853231
    ## [40,] 0.1414214  0.12427894 -0.015197121 -0.14853231
    ## [41,] 0.1414214  0.12427894 -0.015197121 -0.14853231
    ## [42,] 0.1414214  0.12427894 -0.015197121 -0.14853231
    ## [43,] 0.1414214  0.12427894 -0.015197121 -0.14853231
    ## [44,] 0.1414214  0.17831326  0.087358251 -0.06480168
    ## [45,] 0.1414214  0.20533043  0.151676509  0.02419732
    ## [46,] 0.1414214  0.23234759  0.224688481  0.15056334
    ## [47,] 0.1414214  0.23234759  0.224688481  0.15056334
    ## [48,] 0.1414214  0.23234759  0.224688481  0.15056334
    ## [49,] 0.1414214  0.23234759  0.224688481  0.15056334
    ## [50,] 0.1414214  0.25936475  0.306394168  0.31875479

    gammaspeed <- crossprod(gramspeed3$Q,dist)
    gammaspeed

    ##           [,1]
    ## [1,] 303.91449
    ## [2,] 145.55226
    ## [3,]  22.99576
    ## [4,]  13.79688

    coef(lm(dist ~ gramspeed3$Q - 1))

    ## gramspeed3$Q1 gramspeed3$Q2 gramspeed3$Q3 gramspeed3$Q4 
    ##     303.91449     145.55226      22.99576      13.79688

1.  

<!-- -->

    inv.upper.tri(gramspeed3$R,F) %*% gammaspeed

    ##              [,1]
    ## [1,] -19.50504910
    ## [2,]   6.80110598
    ## [3,]  -0.34965781
    ## [4,]   0.01025205

    coef(lm(dist ~ vandermonde(speed,3) - 1))

    ## vandermonde(speed, 3)1 vandermonde(speed, 3)2 vandermonde(speed, 3)3 
    ##           -19.50504910             6.80110597            -0.34965781 
    ## vandermonde(speed, 3)4 
    ##             0.01025205

iii

    summary(lm(dist ~ gramspeed3$Q - 1))

    ## 
    ## Call:
    ## lm(formula = dist ~ gramspeed3$Q - 1)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -26.670  -9.601  -2.231   7.075  44.691 
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)    
    ## gramspeed3$Q1    303.9       15.2  19.988  < 2e-16 ***
    ## gramspeed3$Q2    145.6       15.2   9.573  1.6e-12 ***
    ## gramspeed3$Q3     23.0       15.2   1.512    0.137    
    ## gramspeed3$Q4     13.8       15.2   0.907    0.369    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 15.2 on 46 degrees of freedom
    ## Multiple R-squared:  0.9149, Adjusted R-squared:  0.9075 
    ## F-statistic: 123.6 on 4 and 46 DF,  p-value: < 2.2e-16

    summary(lm(dist ~ vandermonde(speed,3) - 1))

    ## 
    ## Call:
    ## lm(formula = dist ~ vandermonde(speed, 3) - 1)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -26.670  -9.601  -2.231   7.075  44.691 
    ## 
    ## Coefficients:
    ##                         Estimate Std. Error t value Pr(>|t|)
    ## vandermonde(speed, 3)1 -19.50505   28.40530  -0.687    0.496
    ## vandermonde(speed, 3)2   6.80111    6.80113   1.000    0.323
    ## vandermonde(speed, 3)3  -0.34966    0.49988  -0.699    0.488
    ## vandermonde(speed, 3)4   0.01025    0.01130   0.907    0.369
    ## 
    ## Residual standard error: 15.2 on 46 degrees of freedom
    ## Multiple R-squared:  0.9149, Adjusted R-squared:  0.9075 
    ## F-statistic: 123.6 on 4 and 46 DF,  p-value: < 2.2e-16

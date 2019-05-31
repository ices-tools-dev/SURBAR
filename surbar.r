#############################################
# SURBA R implementation: run-time code     #
#                                           #
# Coby Needle, Marine Laboratory            #
#                                           #
# Uses R 3.4.4 (64-bit)	                    # 
#                                           #
# Modifications:                            #
# 25 October 2009 - added fixed parameters, #
# worked on improving speed                 #
# 2 November 2009 - looked again at         #
# dimensions                                #
# 6 November 2009 - coded minimiser from    #
# minpack.lm library                        #
# 26 September 2010 - recoded for NS        #
# haddock example (WKADSAM)                 #
# 7 December 2010 - updated to include      #
# multiple age-structured tuning series     #
# 14 December 2010 - modularised code and   #
# renamed                                   #
# 10 January 2012 - updated version for     #
# WGISDAA, Copenhagen                       #
# 8 February 2012 - updated version for     #
# WKROUND, Aberdeen                         #
# 30 April 2012 - updated version for       #
# WGNSSK, Copenhagen                        #
# 12 March 2014 - updated version for WKHAD #
# 28 April 2017 - updated version for       #
# WGNSSK, Copenhagen                        #
# 5 February 2018 - updated version for     #
# WKNSEA, Copenhagen                        #
# 26 April 2018 - small updates for         #
# WGNSSK, Oostende                          #
# 24 April 2019 - small updates for         #
# WGNSSK, Bergen                            #
# 31 May 2019 - modified for GitHub         #
# repository and TAF                        #
#############################################

# Clear user workspace before beginning
# rm(list=ls(pos=1))

# Close existing graphics windows
# graphics.off()

# Libraries
library(minpack.lm)
library(MASS)
library(grid)
library(lattice)

# Read SURBAR functions
path <- getwd()
source("surbar_functions.r")

# Determine output destination
pdf.plots <- FALSE
x.win <- 10
y.win <- 7
if (pdf.plots)
{
	pdf(file = "surbar_output.pdf",
		onefile = TRUE, width = x.win, height = y.win)
}

# Start timer
start.time <- Sys.time()

# Set up sensitivity parameters (if required)
s.surv <- c("q1q3") # c("q1q3", "q1", "q3")
s.ref <- c(3) # c(2,3)
s.lambda <- c(3) # c(1,3,5)
s.q <- c(1.0) # c(1.0,0.75,0.5)
s.list <- vector("list", length = length(s.surv) * length(s.ref) * length(s.lambda) *
	length(s.q))

# Loop over sensitivity parameters
i.list <- 0
for (i.surv in s.surv)
{
	for (i.ref in s.ref)
	{
		for (i.lambda in s.lambda)
		{
			for (i.q in s.q)
			{
	
i.surv <- s.surv
i.ref <- s.ref
i.lambda <- s.lambda
i.q <- s.q
		
# North Sea lemon sole WGNSSK 2018
sw.file <- "ns_lem_sw.dat"
mat.file <- "ns_lem_mo.dat"
survey.file <- "ns_lem_gam_noage1q1.dat"
zbarage.1 <- 3 
zbarage.2 <- 5 
ref.age <- i.ref 
lambda <- i.lambda
spawn.time <- 0.5
wt.vec <- rep(1.0, 9)

# Read in stock weight data
sw.temp <- read.vpa.file(sw.file)
sw <- data.frame(sw.temp$tab)
sw.a1 <- sw.temp$a[1]
sw.a2 <- sw.temp$a[2]
sw.na <- length(sw.a1:sw.a2)
sw.y1 <- sw.temp$y[1]
sw.y2 <- sw.temp$y[2]
sw.ny <- length(sw.y1:sw.y2)
names(sw) <- sw.a1:sw.a2
rownames(sw) <- sw.y1:sw.y2

# Read in maturity data
mat.temp <- read.vpa.file(mat.file)
mat <- data.frame(mat.temp$tab)
mat.a1 <- mat.temp$a[1]
mat.a2 <- mat.temp$a[2]
mat.na <- length(mat.a1:mat.a2)
mat.y1 <- mat.temp$y[1]
mat.y2 <- mat.temp$y[2]
mat.ny <- length(mat.y1:mat.y2)
names(mat) <- mat.a1:mat.a2
rownames(mat) <- mat.y1:mat.y2

# Read in survey data
idx.temp <- read.survey.file(survey.file)
numk <- idx.temp$n
idx <- idx.temp$idx
rho <- unlist(lapply(idx, function(wk){wk$rho}))
cat("\nIndices read:")
cat(paste("\n", unlist(lapply(idx, function(wk){wk$name})), sep = ""))
cat("\n")

##################################################
# Catchability adjustments should be applied here!
##################################################

# Trim survey dataset according to sensitivity parameter
if (i.surv == "q1")
{
	numk <- 1
	idx <- list(idx[[1]])
	rho <- rho[1]
} else if (i.surv == "q3")
{
	numk <- 1
	idx <- list(idx[[2]])
	rho <- rho[2]
}

cat("\nIndices used:")
cat(paste("\n", unlist(lapply(idx, function(wk){wk$name})), sep = ""))
cat("\n")

cat("\nSettings:")
cat(paste0("\nReference age = ", i.ref))
cat(paste0("\nLambda = ", i.lambda))
cat(paste0("\nCatchability = ", i.q))

# Set survey data to consistent dimensions
x <- vector("list", length = numk)
names(x) <- unlist(lapply(idx, function(wk){wk$name}))
x <- lapply(x, function(wk)
	{
		wk.tab <- data.frame(array(NA, dim = dim(sw)))
		names(wk.tab) <- names(sw)
		rownames(wk.tab) <- rownames(sw)
		wk.tab
	})
for (k in 1:numk)
{
	a <- idx[[k]]$a1:idx[[k]]$a2
	y <- idx[[k]]$y1:idx[[k]]$y2
	x[[k]][ch(y),ch(a)] <- idx[[k]]$tab[ch(y),ch(a)]
}

# Replace zeros with minimum value
for (k in 1:numk)
{
	xx <- x[[k]]
	x.min1 <- min(xx, na.rm = TRUE)
	if (x.min1 == 0)
	{
		xx[xx == 0] <- 9999
		x.min2 <- min(xx, na.rm = TRUE)
		xx[xx == 9999.0] <- x.min2
	}
	x[[k]] <- xx
}

# Mean-standardise survey data
for (k in 1:numk)
{
	x.mean <- mean(apply(x[[k]], 2, function(wk){mean(wk, na.rm = TRUE)}), na.rm = TRUE)
	x[[k]] <- x[[k]] / x.mean
}

# Set up catchability and weightings arrays
# 1. Reduce catchability on youngest ages (hook in catch curves)
# 2. Reduce weighting on points with |residuals| > 2
qval <- vector("list", length = numk)
wt <- vector("list", length = numk)
for (k in 1:numk)
{
	qval[[k]] <- data.frame(array(1.0, dim = dim(x[[k]])))
	rownames(qval[[k]]) <- rownames(x[[k]])
	names(qval[[k]]) <- names(x[[k]])

	qval[[k]][,1] <- 0.1
	qval[[k]][,2] <- 0.5
	qval[[k]][,6:9] <- i.q
	
	wt[[k]] <- data.frame(array(wt.vec[k], dim = dim(x[[k]])))
	rownames(wt[[k]]) <- rownames(x[[k]])
	names(wt[[k]]) <- names(x[[k]])
}

# Alternatives:
# Q1 catchability estimates (from spreadsheet analysis)
#qval[[1]][,1] <- 0.001
#qval[[1]][,2] <- 0.063802325
#qval[[1]][,3] <- 0.360767917
#qval[[1]][,6:9] <- i.q

#qval[[1]][,1] <- 0.01
#qval[[1]][,2] <- 0.5
#qval[[1]][,6:9] <- i.q

# Q3 catchability estimates (from spreadsheet analysis)
#qval[[2]][,1] <- 0.001
#qval[[2]][,2] <- 0.059161712
#qval[[2]][,3] <- 0.161156136
#qval[[2]][,6:9] <- i.q

#qval[[2]][,1] <- 0.01
#qval[[2]][,2] <- 0.5
#qval[[2]][,6:9] <- i.q

# Optional: downweighting select points
#wt[[1]][as.character(2008),1] <- 0.5
#wt[[2]][as.character(2013),1] <- 0.5

# Set up maximal age and year ranges across surveys and trim accordingly
x.temp <- trim.data(x, idx, qval, wt, sw, mat)
x <- x.temp$x
qval <- x.temp$qval
wt <- x.temp$wt
sw <- x.temp$sw
mat <- x.temp$mat

a1 <- min(as.numeric(names(sw)))
a2 <- max(as.numeric(names(sw)))
na <- length(a1:a2)
y1 <- min(as.numeric(rownames(sw)))
y2 <- max(as.numeric(rownames(sw)))
ny <- length(y1:y2)

# General target vectors (leaving room for fixed values)

s0.a <- seq(1.0, 1.0, length = ref.age - a1 + 1)
s0.b <- seq(1.0, 1.0, length = a2 - ref.age)
s0 <- c(s0.a[1:(length(s0.a) - 1)], s0.b[2:length(s0.b)])
if (length(na.omit(s0)) < length(s0))
{
	stop("\nCheck reference age.\n")
}
f0 <- rep(1.0, ny - 1)

# Initial estimates for r are taken from the log of the average (across surveys)
# of the survey index values at the appropriate years/ages
temp.r <- data.frame(array(NA, dim = dim(sw)))
names(temp.r) <- names(sw)
rownames(temp.r) <- rownames(sw)
for (i in rownames(temp.r))
{
	for (j in names(temp.r))
	{
		temp.r[i,j] <- mean(unlist(lapply(x, function(wk){wk[i, j]})), na.rm = TRUE)
	}
}
# Fill missing values with age-based averages
for (j in names(temp.r))
{
	temp.mean <- mean(temp.r[,j], na.rm = TRUE)
	temp.r[is.nan(temp.r[,j]),j] <- temp.mean
}
r0 <- log(as.numeric(unlist(c(rev(temp.r[1,]), temp.r[2:dim(temp.r)[1],1]))))
params0 <- c(s0, f0, r0)

# Test: if there are ages and years for which there are NO data in the initial estimates of r, stop
if (length(is.nan(r0)[is.nan(r0)]) > 0)
{
	cat("\nAcross-survey averages:\n")
	print(temp.r)
	stop("One of the recruiting cohorts lacks information from any survey.")
}

# Set up control object

mqdt.control <- list(ftol = 0.00001, ptol = 0.00001, gtol = 0, 
	diag = numeric(), factor = 100, maxfev = 100 * (length(params0) + 1), 
	nprint = 1, maxiter = 200)

# Set bounds on parameters
s0.lower <- rep(-5, length = length(s0))
s0.upper <- rep(5, length = length(s0))
f0.lower <- rep(-1, length = length(f0))
f0.upper <- rep(5, length = length(f0))
r0.lower <- rep(-Inf, length = length(r0))
r0.upper <- rep(Inf, length = length(r0))
lower0 <- c(s0.lower, f0.lower, r0.lower)
upper0 <- c(s0.upper, f0.upper, r0.upper)

# Run minimisation
cat("\nOptimising parameters...\n")
x.surba <- nls.lm(params0, lower = lower0, upper = upper0, 
	fn = surbar, control = mqdt.control)

# Extract parameter values

s <- rep(NA, length = na)
s[1:(ref.age-1)] <- x.surba$par[1:(ref.age-1)]
s[ref.age] <- 1.0
s[(ref.age+1):(na-1)] <- x.surba$par[ref.age:(na-2)]
s[na] <- s[na-1]

f <- rep(NA, length = ny)
f[1:(ny-1)] <- x.surba$par[(na-1):(na + ny - 3)]
f[ny] <- mean(f[(ny-3):(ny-1)])

r <- x.surba$par[(na + ny - 2):length(x.surba$par)]

# Extract log residuals

res <- vector("list", length = numk)
x.res <- surbar(x.surba$par)
res.start <- 1
for (k in 1:numk)
{
	x.temp <- x[[k]][rownames(idx[[k]]$tab),names(idx[[k]]$tab)]
	res.end <- (res.start - 1) + (dim(x.temp)[1] * dim(x.temp)[2])
	res[[k]] <- data.frame(array(x.res[res.start:res.end], dim = dim(x.temp)))
	names(res[[k]]) <- names(x.temp)
	rownames(res[[k]]) <- rownames(x.temp)
	res.start <- res.end + 1
}

# Construct population estimates

stock <- data.frame(array(NA, dim = c(ny,5)))
names(stock) <- c("year", "rec", "ssb", "tsb", "meanz")

zmort <- f %o% s

lnn <- array(NA, dim = dim(zmort))
lnn[1,] <- rev(r[1:dim(lnn)[2]])
lnn[2:dim(lnn)[1],1] <- r[(dim(lnn)[2]+1):length(r)]
# Oldest true age
wk.ota <- dim(zmort)[1]
for (j in 2:dim(zmort)[2])
{
	for (i in 2:dim(zmort)[1])
	{
		lnn[i,j] <- lnn[i-1,j-1] - zmort[i-1,j-1]
	}
}
n <- exp(lnn)	

stock$year <- y1:y2
stock$meanz <- apply(zmort[,zbarage.1:zbarage.2], 1, mean)
stock$rec <- exp(r[na:length(r)])
stock$ssb <- apply(n * sw * mat, 1, sum)	
stock$tsb <- apply(n * sw, 1, sum)

# Mean-standardise rec, ssb, tsb (if required)
#mean.rec <- mean(stock$rec)
#mean.ssb <- mean(stock$ssb)
#mean.tsb <- mean(stock$tsb)
#stock$rec <- stock$rec / mean.rec
#stock$ssb <- stock$ssb / mean.ssb
#stock$tsb <- stock$tsb / mean.tsb

cat("\nStock summary:\n")
print(round(stock,3))

i.list <- i.list + 1
s.list[[i.list]] <- list(surv = i.surv, ref = i.ref, lambda = i.lambda,
	q = i.q, dev = deviance(x.surba), aic = AIC.nls.lm(x.surba), stock = stock)
cat(paste0("Completed ", i.list, " of ", length(s.list)))

}}}} # End sensitivity loop

# Sensitivity analysis

if(FALSE){

# Effect of survey choice
list.condition <- sapply(s.list, function(x){x$ref == 3 & x$lambda == 3.0 & x$q == 1})
s0.list  <- s.list[list.condition]
if (!pdf.plots) windows(width = 10, height = 7)
par(mfrow = c(2,2), mar = c(3,4,2,1))
# Mean Z
ymin <- min(0,unlist(lapply(s0.list, function(wk){min(wk$stock$meanz)})))
ymax <- max(unlist(lapply(s0.list, function(wk){max(wk$stock$meanz)})))
plot(s0.list[[1]]$stock$year, s0.list[[1]]$stock$meanz, type = "l", lty = 1, col = 1,
	xlab = "", ylab = "Mean Z", ylim = c(ymin,ymax))
for (i in 2:length(s0.list))
{
	lines(s0.list[[i]]$stock$year, s0.list[[i]]$stock$meanz, type = "l", lty = 1, col = i)
}
abline(h = 0, lty = 8)
legend(x = "topright", bty = "n", legend = s.surv, lty = 1, col = 1:length(s0.list))
# SSB
ymin <- min(0,unlist(lapply(s0.list, function(wk){min(wk$stock$ssb)})))
ymax <- max(unlist(lapply(s0.list, function(wk){max(wk$stock$ssb)})))
plot(s0.list[[1]]$stock$year, s0.list[[1]]$stock$ssb, type = "l", lty = 1, col = 1,
	xlab = "", ylab = "SSB", ylim = c(ymin,ymax))
for (i in 2:length(s0.list))
{
	lines(s0.list[[i]]$stock$year, s0.list[[i]]$stock$ssb, type = "l", lty = 1, col = i)
}
abline(h = 0, lty = 8)
legend(x = "topright", bty = "n", legend = s.surv, lty = 1, col = 1:length(s0.list))
# TSB
ymin <- min(0,unlist(lapply(s0.list, function(wk){min(wk$stock$tsb)})))
ymax <- max(unlist(lapply(s0.list, function(wk){max(wk$stock$tsb)})))
plot(s0.list[[1]]$stock$year, s0.list[[1]]$stock$tsb, type = "l", lty = 1, col = 1,
	xlab = "", ylab = "TSB", ylim = c(ymin,ymax))
for (i in 2:length(s0.list))
{
	lines(s0.list[[i]]$stock$year, s0.list[[i]]$stock$tsb, type = "l", lty = 1, col = i)
}
abline(h = 0, lty = 8)
legend(x = "topright", bty = "n", legend = s.surv, lty = 1, col = 1:length(s0.list))
# Recruitment
ymin <- min(0,unlist(lapply(s0.list, function(wk){min(wk$stock$rec)})))
ymax <- max(unlist(lapply(s0.list, function(wk){max(wk$stock$rec)})))
plot(s0.list[[1]]$stock$year, s0.list[[1]]$stock$rec, type = "l", lty = 1, col = 1,
	xlab = "", ylab = "Rec at age 1", ylim = c(ymin,ymax))
for (i in 2:length(s0.list))
{
	lines(s0.list[[i]]$stock$year, s0.list[[i]]$stock$rec, type = "l", lty = 1, col = i)
}
abline(h = 0, lty = 8)
legend(x = "topright", bty = "n", legend = s.surv, lty = 1, col = 1:length(s0.list))
mtext("Effect of survey used", side = 3, outer = TRUE, line = -1.5)

# Effect of reference age
list.condition <- sapply(s.list, function(x){x$surv == "q1q3" & x$lambda == 3.0 & x$q == 1})
s0.list  <- s.list[list.condition]
if (!pdf.plots) windows(width = 10, height = 7)
par(mfrow = c(2,2), mar = c(3,4,2,1))
# Mean Z
ymin <- min(0,unlist(lapply(s0.list, function(wk){min(wk$stock$meanz)})))
ymax <- max(unlist(lapply(s0.list, function(wk){max(wk$stock$meanz)})))
plot(s0.list[[1]]$stock$year, s0.list[[1]]$stock$meanz, type = "l", lty = 1, col = 1,
	xlab = "", ylab = "Mean Z", ylim = c(ymin,ymax))
for (i in 2:length(s0.list))
{
	lines(s0.list[[i]]$stock$year, s0.list[[i]]$stock$meanz, type = "l", lty = 1, col = i)
}
abline(h = 0, lty = 8)
legend(x = "topright", bty = "n", legend = s.ref, lty = 1, col = 1:length(s0.list))
# SSB
ymin <- min(0,unlist(lapply(s0.list, function(wk){min(wk$stock$ssb)})))
ymax <- max(unlist(lapply(s0.list, function(wk){max(wk$stock$ssb)})))
plot(s0.list[[1]]$stock$year, s0.list[[1]]$stock$ssb, type = "l", lty = 1, col = 1,
	xlab = "", ylab = "SSB", ylim = c(ymin,ymax))
for (i in 2:length(s0.list))
{
	lines(s0.list[[i]]$stock$year, s0.list[[i]]$stock$ssb, type = "l", lty = 1, col = i)
}
abline(h = 0, lty = 8)
legend(x = "topright", bty = "n", legend = s.ref, lty = 1, col = 1:length(s0.list))
# TSB
ymin <- min(0,unlist(lapply(s0.list, function(wk){min(wk$stock$tsb)})))
ymax <- max(unlist(lapply(s0.list, function(wk){max(wk$stock$tsb)})))
plot(s0.list[[1]]$stock$year, s0.list[[1]]$stock$tsb, type = "l", lty = 1, col = 1,
	xlab = "", ylab = "TSB", ylim = c(ymin,ymax))
for (i in 2:length(s0.list))
{
	lines(s0.list[[i]]$stock$year, s0.list[[i]]$stock$tsb, type = "l", lty = 1, col = i)
}
abline(h = 0, lty = 8)
legend(x = "topright", bty = "n", legend = s.ref, lty = 1, col = 1:length(s0.list))
# Recruitment
ymin <- min(0,unlist(lapply(s0.list, function(wk){min(wk$stock$rec)})))
ymax <- max(unlist(lapply(s0.list, function(wk){max(wk$stock$rec)})))
plot(s0.list[[1]]$stock$year, s0.list[[1]]$stock$rec, type = "l", lty = 1, col = 1,
	xlab = "", ylab = "Rec at age 1", ylim = c(ymin,ymax))
for (i in 2:length(s0.list))
{
	lines(s0.list[[i]]$stock$year, s0.list[[i]]$stock$rec, type = "l", lty = 1, col = i)
}
abline(h = 0, lty = 8)
legend(x = "topright", bty = "n", legend = s.ref, lty = 1, col = 1:length(s0.list))
mtext("Effect of reference age", side = 3, outer = TRUE, line = -1.5)

# Effect of lambda
list.condition <- sapply(s.list, function(x){x$surv == "q1q3" & x$ref == 3 & x$q == 1})
s0.list  <- s.list[list.condition]
if (!pdf.plots) windows(width = 10, height = 7)
par(mfrow = c(2,2), mar = c(3,4,2,1))
# Mean Z
ymin <- min(0,unlist(lapply(s0.list, function(wk){min(wk$stock$meanz)})))
ymax <- max(unlist(lapply(s0.list, function(wk){max(wk$stock$meanz)})))
plot(s0.list[[1]]$stock$year, s0.list[[1]]$stock$meanz, type = "l", lty = 1, col = 1,
	xlab = "", ylab = "Mean Z", ylim = c(ymin,ymax))
for (i in 2:length(s0.list))
{
	lines(s0.list[[i]]$stock$year, s0.list[[i]]$stock$meanz, type = "l", lty = 1, col = i)
}
abline(h = 0, lty = 8)
legend(x = "topright", bty = "n", legend = s.lambda, lty = 1, col = 1:length(s0.list))
# SSB
ymin <- min(0,unlist(lapply(s0.list, function(wk){min(wk$stock$ssb)})))
ymax <- max(unlist(lapply(s0.list, function(wk){max(wk$stock$ssb)})))
plot(s0.list[[1]]$stock$year, s0.list[[1]]$stock$ssb, type = "l", lty = 1, col = 1,
	xlab = "", ylab = "SSB", ylim = c(ymin,ymax))
for (i in 2:length(s0.list))
{
	lines(s0.list[[i]]$stock$year, s0.list[[i]]$stock$ssb, type = "l", lty = 1, col = i)
}
abline(h = 0, lty = 8)
legend(x = "topright", bty = "n", legend = s.lambda, lty = 1, col = 1:length(s0.list))
# TSB
ymin <- min(0,unlist(lapply(s0.list, function(wk){min(wk$stock$tsb)})))
ymax <- max(unlist(lapply(s0.list, function(wk){max(wk$stock$tsb)})))
plot(s0.list[[1]]$stock$year, s0.list[[1]]$stock$tsb, type = "l", lty = 1, col = 1,
	xlab = "", ylab = "TSB", ylim = c(ymin,ymax))
for (i in 2:length(s0.list))
{
	lines(s0.list[[i]]$stock$year, s0.list[[i]]$stock$tsb, type = "l", lty = 1, col = i)
}
abline(h = 0, lty = 8)
legend(x = "topright", bty = "n", legend = s.lambda, lty = 1, col = 1:length(s0.list))
# Recruitment
ymin <- min(0,unlist(lapply(s0.list, function(wk){min(wk$stock$rec)})))
ymax <- max(unlist(lapply(s0.list, function(wk){max(wk$stock$rec)})))
plot(s0.list[[1]]$stock$year, s0.list[[1]]$stock$rec, type = "l", lty = 1, col = 1,
	xlab = "", ylab = "Rec at age 1", ylim = c(ymin,ymax))
for (i in 2:length(s0.list))
{
	lines(s0.list[[i]]$stock$year, s0.list[[i]]$stock$rec, type = "l", lty = 1, col = i)
}
abline(h = 0, lty = 8)
legend(x = "topright", bty = "n", legend = s.lambda, lty = 1, col = 1:length(s0.list))
mtext("Effect of smoother", side = 3, outer = TRUE, line = -1.5)

# Effect of q
list.condition <- sapply(s.list, function(x){x$surv == "q1q3" & x$ref == 3 & x$lambda == 3})
s0.list  <- s.list[list.condition]
if (!pdf.plots) windows(width = 10, height = 7)
par(mfrow = c(2,2), mar = c(3,4,2,1))
# Mean Z
ymin <- min(0,unlist(lapply(s0.list, function(wk){min(wk$stock$meanz)})))
ymax <- max(unlist(lapply(s0.list, function(wk){max(wk$stock$meanz)})))
plot(s0.list[[1]]$stock$year, s0.list[[1]]$stock$meanz, type = "l", lty = 1, col = 1,
	xlab = "", ylab = "Mean Z", ylim = c(ymin,ymax))
for (i in 2:length(s0.list))
{
	lines(s0.list[[i]]$stock$year, s0.list[[i]]$stock$meanz, type = "l", lty = 1, col = i)
}
abline(h = 0, lty = 8)
legend(x = "topright", bty = "n", legend = s.q, lty = 1, col = 1:length(s0.list))
# SSB
ymin <- min(0,unlist(lapply(s0.list, function(wk){min(wk$stock$ssb)})))
ymax <- max(unlist(lapply(s0.list, function(wk){max(wk$stock$ssb)})))
plot(s0.list[[1]]$stock$year, s0.list[[1]]$stock$ssb, type = "l", lty = 1, col = 1,
	xlab = "", ylab = "SSB", ylim = c(ymin,ymax))
for (i in 2:length(s0.list))
{
	lines(s0.list[[i]]$stock$year, s0.list[[i]]$stock$ssb, type = "l", lty = 1, col = i)
}
abline(h = 0, lty = 8)
legend(x = "topright", bty = "n", legend = s.q, lty = 1, col = 1:length(s0.list))
# TSB
ymin <- min(0,unlist(lapply(s0.list, function(wk){min(wk$stock$tsb)})))
ymax <- max(unlist(lapply(s0.list, function(wk){max(wk$stock$tsb)})))
plot(s0.list[[1]]$stock$year, s0.list[[1]]$stock$tsb, type = "l", lty = 1, col = 1,
	xlab = "", ylab = "TSB", ylim = c(ymin,ymax))
for (i in 2:length(s0.list))
{
	lines(s0.list[[i]]$stock$year, s0.list[[i]]$stock$tsb, type = "l", lty = 1, col = i)
}
abline(h = 0, lty = 8)
legend(x = "topright", bty = "n", legend = s.q, lty = 1, col = 1:length(s0.list))
# Recruitment
ymin <- min(0,unlist(lapply(s0.list, function(wk){min(wk$stock$rec)})))
ymax <- max(unlist(lapply(s0.list, function(wk){max(wk$stock$rec)})))
plot(s0.list[[1]]$stock$year, s0.list[[1]]$stock$rec, type = "l", lty = 1, col = 1,
	xlab = "", ylab = "Rec at age 1", ylim = c(ymin,ymax))
for (i in 2:length(s0.list))
{
	lines(s0.list[[i]]$stock$year, s0.list[[i]]$stock$rec, type = "l", lty = 1, col = i)
}
abline(h = 0, lty = 8)
legend(x = "topright", bty = "n", legend = s.q, lty = 1, col = 1:length(s0.list))
mtext("Effect of catchability", side = 3, outer = TRUE, line = -1.5)

#if (pdf.plots) dev.off()
#stop.time <- Sys.time()
#cat(paste("\nRun complete in:", round(stop.time - start.time, 3), "s\n"))
#stop() # Stop here when only doing sensitivity loop

} # End sensitivity analysis output

#####################################################################################

# Bootstrap model fit using data simulation

cat("\nBootstrapping parameters...\n")

# Test for singularity
x.eigen <- eigen(x.surba$hessian, only.values = TRUE)$values
if (length(x.eigen[x.eigen == 0]) > 0)
{
	stop("At least one parameter cannot be estimated: \nCheck that there are data for each cohort.")
}

n.psim <- 1000
x.psim <- mvrnorm(n = n.psim, mu = x.surba$par, vcov(x.surba))

x.stock <- data.frame(array(NA, dim = c(ny, 5)))
names(x.stock) <- c("year", "rec", "ssb", "tsb", "meanz")
x.psim.stock <- vector("list", n.psim)
x.psim.stock <- lapply(x.psim.stock, function(wk){wk <- x.stock})

x.psim.s <- array(NA, dim = c(1000, na))
x.psim.s[,1:(ref.age-1)] <- x.psim[,1:(ref.age-1)]
x.psim.s[,ref.age] <- 1.0
x.psim.s[,(ref.age+1):(na-1)] <- x.psim[,ref.age:(na-2)]
x.psim.s[,na] <- x.psim.s[,na-1]

x.psim.f <- array(NA, dim = c(1000, ny))
x.psim.f[,1:(ny-1)] <- x.psim[,(na-1):(na + ny - 3)]
x.psim.f[,ny] <- apply(x.psim.f[,(ny-3):(ny-1)], 1, mean)

x.psim.r <- x.psim[,(na + ny - 2):length(x.surba$par)]

x.s <- rep(NA, length = na)
x.f <- rep(NA, length = ny)

for (i in 1:n.psim)
{
	# Extract parameters

	x.s[1:(ref.age-1)] <- x.psim[i,1:(ref.age-1)]
	x.s[ref.age] <- 1.0
	x.s[(ref.age+1):(na-1)] <- x.psim[i,ref.age:(na-2)]
	x.s[na] <- x.s[na-1]

	x.f[1:(ny-1)] <- x.psim[i,(na-1):(na + ny - 3)]
	x.f[ny] <- mean(x.f[(ny-3):(ny-1)])

	x.r <- x.psim[i,(na + ny - 2):length(x.surba$par)]

	# Generate stock summaries

	zmort <- x.f %o% x.s

	lnn <- array(NA, dim = dim(zmort))
	lnn[1,] <- rev(x.r[1:dim(lnn)[2]])
	lnn[2:dim(lnn)[1],1] <- x.r[(dim(lnn)[2]+1):length(x.r)]
	for (jj in 2:dim(zmort)[2])
	{
		for (ii in 2:dim(zmort)[1])
		{
			lnn[ii,jj] <- lnn[ii-1,jj-1] - zmort[ii-1,jj-1]
		}
	}
	n <- exp(lnn)	

	x.psim.stock[[i]]$year <- y1:y2
	x.psim.stock[[i]]$meanz <- apply(zmort[,zbarage.1:zbarage.2], 1, mean)
	x.psim.stock[[i]]$z <- zmort
	x.psim.stock[[i]]$rec <- exp(x.r[na:length(r)]) # / mean.rec
	x.psim.stock[[i]]$ssb <- apply(n * sw * mat, 1, sum) # / mean.ssb
	x.psim.stock[[i]]$tsb <- apply(n * sw, 1, sum) # / mean.tsb
}

# Tabulate lower and upper bounds of estimates (mean Z, SSB, recruitment)

# Mean Z
x.psim.meanz <- do.call(rbind, lapply(x.psim.stock, function(wk){wk$meanz}))
x.psim.meanz.quantile <- array(NA, dim = c(dim(x.psim.meanz)[2],5))
x.psim.meanz.mean <- rep(NA, dim(x.psim.meanz)[2])
rownames(x.psim.meanz.quantile) <- y1:y2
colnames(x.psim.meanz.quantile) <- c("5%","25%","50%","75%","95%")
for (i in 1:dim(x.psim.meanz)[2])
{
	x.psim.meanz.quantile[i,] <- quantile(x.psim.meanz[,i], c(0.05, 0.25, 0.5, 0.75, 0.95))
      x.psim.meanz.mean[i] <- mean(x.psim.meanz[,i])
}

# SSB
x.psim.ssb <- do.call(rbind, lapply(x.psim.stock, function(wk){wk$ssb}))
x.psim.ssb.quantile <- array(NA, dim = c(dim(x.psim.ssb)[2],5))
x.psim.ssb.mean <- rep(NA, dim(x.psim.ssb)[2])
rownames(x.psim.ssb.quantile) <- y1:y2
colnames(x.psim.ssb.quantile) <- c("5%","25%","50%","75%","95%")
for (i in 1:dim(x.psim.ssb)[2])
{
	x.psim.ssb.quantile[i,] <- quantile(x.psim.ssb[,i], c(0.05, 0.25, 0.5, 0.75, 0.95))
      x.psim.ssb.mean[i] <- mean(x.psim.ssb[,i])
}

# Recruitment
x.psim.rec <- do.call(rbind, lapply(x.psim.stock, function(wk){wk$rec}))
x.psim.rec.quantile <- array(NA, dim = c(dim(x.psim.rec)[2],5))
x.psim.rec.mean <- rep(NA, dim(x.psim.rec)[2])
rownames(x.psim.rec.quantile) <- y1:y2
colnames(x.psim.rec.quantile) <- c("5%","25%","50%","75%","95%")
for (i in 1:dim(x.psim.rec)[2])
{
	x.psim.rec.quantile[i,] <- quantile(x.psim.rec[,i], c(0.05, 0.25, 0.5, 0.75, 0.95))
      x.psim.rec.mean[i] <- mean(x.psim.rec[,i])
}

x.psim.summary <- data.frame(array(NA, dim = c(length(y1:y2), 9)))
rownames(x.psim.summary) <- y1:y2
colnames(x.psim.summary) <- c("meanz.lower", "meanz", "meanz.upper", "ssb.lower", 
	"ssb", "ssb.upper", "rec.lower", "rec", "rec.upper")
x.psim.summary$meanz.lower <- as.numeric(x.psim.meanz.quantile[,"5%"])
x.psim.summary$meanz <- stock$meanz
x.psim.summary$meanz.upper <- as.numeric(x.psim.meanz.quantile[,"95%"])
x.psim.summary$ssb.lower <- as.numeric(x.psim.ssb.quantile[,"5%"])
x.psim.summary$ssb <- stock$ssb
x.psim.summary$ssb.upper <- as.numeric(x.psim.ssb.quantile[,"95%"])
x.psim.summary$rec.lower <- as.numeric(x.psim.rec.quantile[,"5%"])
x.psim.summary$rec <- stock$rec
x.psim.summary$rec.upper <- as.numeric(x.psim.rec.quantile[,"95%"])

cat("\nExtended stock summary:\n")
print(round(x.psim.summary, 3))

# Estimate standard errors and CVs of estimates

x.se <- data.frame(array(NA, dim = c(ny, 5)))
names(x.se) <- c("year","meanz","rec","ssb","tsb")
x.se$year <- y1:y2
x.se$meanz <- apply(do.call(rbind, lapply(x.psim.stock, function(wk){wk$meanz})), 2, sd)
x.se$rec <- apply(do.call(rbind, lapply(x.psim.stock, function(wk){wk$rec})), 2, sd)
x.se$ssb <- apply(do.call(rbind, lapply(x.psim.stock, function(wk){wk$ssb})), 2, sd)
x.se$tsb <- apply(do.call(rbind, lapply(x.psim.stock, function(wk){wk$tsb})), 2, sd)

x.cv <- data.frame(array(NA, dim = c(ny, 5)))
names(x.cv) <- c("year","meanz","rec","ssb","tsb")
x.cv$year <- y1:y2
x.cv$meanz <- x.se$meanz / stock$meanz
x.cv$rec <- x.se$rec / stock$rec
x.cv$ssb <- x.se$ssb / stock$ssb
x.cv$tsb <- x.se$tsb / stock$tsb

#####################################################################################

# Summary plots

plot.surbar(x = x, stock = stock, psim = x.psim.stock, psim.s = x.psim.s, psim.f = x.psim.f, psim.r = x.psim.r, 
	res = res, y1 = y1, y2 = y2, a1 = a1, a2 = a2, type = "sum.line", pdf.plots)
plot.surbar(x = x, stock = stock, psim = x.psim.stock, psim.s = x.psim.s, psim.f = x.psim.f, psim.r = x.psim.r, 
	res = res, y1 = y1, y2 = y2, a1 = a1, a2 = a2, type = "sum.boxplot", pdf.plots)
plot.surbar(x = x, stock = stock, psim = x.psim.stock, psim.s = x.psim.s, psim.f = x.psim.f, psim.r = x.psim.r, 
	res = res, y1 = y1, y2 = y2, a1 = a1, a2 = a2, type = "res.line", pdf.plots)
plot.surbar(x = x, stock = stock, psim = x.psim.stock, psim.s = x.psim.s, psim.f = x.psim.f, psim.r = x.psim.r, 
	res = res, y1 = y1, y2 = y2, a1 = a1, a2 = a2, type = "res.smooth", pdf.plots)
plot.surbar(x = x, stock = stock, psim = x.psim.stock, psim.s = x.psim.s, psim.f = x.psim.f, psim.r = x.psim.r, 
	res = res, y1 = y1, y2 = y2, a1 = a1, a2 = a2, type = "params", pdf.plots)

plot.surbar(x = x, stock = stock, psim = x.psim.stock, psim.s = x.psim.s, psim.f = x.psim.f, psim.r = x.psim.r, 
	res = res, y1 = y1, y2 = y2, a1 = a1, a2 = a2, type = "catch.curve", pdf.plots)
plot.surbar(x = x, stock = stock, psim = x.psim.stock, psim.s = x.psim.s, psim.f = x.psim.f, psim.r = x.psim.r, 
	res = res, y1 = y1, y2 = y2, a1 = a1, a2 = a2, type = "log.by.cohort", pdf.plots)

#plot.surbar(x = x, stock = stock, psim = x.psim.stock, psim.s = x.psim.s, psim.f = x.psim.f, psim.r = x.psim.r, 
#	res = res, y1 = y1, y2 = y2, a1 = a1, a2 = a2, type = "age.scatterplot", pdf.plots)
if (!pdf.plots) windows(width = 10, height = 7)
plot.index.corr(wk.object = list(x[[1]][as.character(2007:2018),2:5]), wk.type = "SURBAR") # Q1
if (!pdf.plots) windows(width = 10, height = 7)
plot.index.corr(wk.object = list(x[[2]][as.character(2006:2018),1:9]), wk.type = "SURBAR") # Q3

#####################################################################################

# Retrospective runs

xr <- x; qvalr <- qval
wtr <- wt; swr <- sw; matr <- mat
y2 <- max(as.numeric(rownames(sw)))
nretro <- 5

retro.stock <- vector("list", length = nretro)
kk <- 0

for (ry in (y2-1):(y2-nretro))
{
	kk <- kk + 1

	# Trim datasets	
	for (k in 1:length(xr))
	{
		xr[[k]] <- x[[k]][row.names(x[[k]]) <= ry,]
	}
	xr.temp <- trim.retro.data(xr, qvalr, wtr, swr, matr)
	xr <- xr.temp$x
	qvalr <- xr.temp$qval
	wtr <- xr.temp$wt
	swr <- xr.temp$sw
	matr <- xr.temp$mat

	# Set up year range for retro run
	y1 <- min(as.numeric(rownames(swr)))
	y2 <- max(as.numeric(rownames(swr)))
	ny <- length(y1:y2)

	# General target vectors (leaving room for fixed values)
	s0.a <- seq(1.0, 1.0, length = ref.age - a1 + 1)
	s0.b <- seq(1.0, 1.0, length = a2 - ref.age)
	s0 <- c(s0.a[1:(length(s0.a) - 1)], s0.b[2:length(s0.b)])
	if (length(na.omit(s0)) < length(s0))
	{
		stop("\nCheck reference age.\n")
	}
	f0 <- rep(1.0, ny - 1)

	# Initial estimates for r are taken from the log of the average (across surveys)
	# of the survey index values at the appropriate years/ages
	temp.r <- data.frame(array(NA, dim = dim(swr)))
	names(temp.r) <- names(swr)
	rownames(temp.r) <- rownames(swr)
	for (i in rownames(temp.r))
	{
		for (j in names(temp.r))
		{
			temp.r[i,j] <- mean(unlist(lapply(xr, function(wk){wk[i, j]})), na.rm = TRUE)
		}
	}
	# Fill missing values with age-based averages
	for (j in names(temp.r))
	{
		temp.mean <- mean(temp.r[,j], na.rm = TRUE)
		temp.r[is.nan(temp.r[,j]),j] <- temp.mean
	}
	r0 <- log(as.numeric(unlist(c(rev(temp.r[1,]), temp.r[2:dim(temp.r)[1],1]))))
	params0 <- c(s0, f0, r0)

	# Test: if there are ages and years for which there are NO data in the initial estimates of r, stop
	if (length(is.nan(r0)[is.nan(r0)]) > 0)
	{
		cat("\nAcross-survey averages:\n")
		print(temp.r)
		stop("One of the recruiting cohorts lacks information from any survey.")
	}

	# Set up control object
	mqdt.control <- list(ftol = 0.00001, ptol = 0.00001, gtol = 0, 
		diag = numeric(), factor = 100, maxfev = 100 * (length(params0) + 1), 
		nprint = 1, maxiter = 200)

	# Run minimisation
	cat(paste("\nOptimising retro parameters - last year ", ry, "...\n", sep = ""))
	xr.surba <- nls.lm(params0, fn = surbar.retro, control = mqdt.control)

	# Extract parameter values
	sr <- rep(NA, length = na)
	sr[1:(ref.age-1)] <- xr.surba$par[1:(ref.age-1)]
	sr[ref.age] <- 1.0
	sr[(ref.age+1):(na-1)] <- xr.surba$par[ref.age:(na-2)]
	sr[na] <- sr[na-1]
	fr <- rep(NA, length = ny)
	fr[1:(ny-1)] <- xr.surba$par[(na-1):(na + ny - 3)]
	fr[ny] <- mean(fr[(ny-3):(ny-1)])
	rr <- xr.surba$par[(na + ny - 2):length(xr.surba$par)]

	# Construct population estimates
	stockr <- data.frame(array(NA, dim = c(ny,5)))
	names(stockr) <- c("year", "rec", "ssb", "tsb", "meanz")
	zmortr <- fr %o% sr
	lnnr <- array(NA, dim = dim(zmortr))
	lnnr[1,] <- rev(rr[1:dim(lnnr)[2]])
	lnnr[2:dim(lnnr)[1],1] <- rr[(dim(lnnr)[2]+1):length(rr)]
	# Oldest true age
	wk.ota <- dim(zmortr)[1]
	for (j in 2:dim(zmortr)[2])
	{
		for (i in 2:dim(zmortr)[1])
		{
			lnnr[i,j] <- lnnr[i-1,j-1] - zmortr[i-1,j-1]
		}
	}
	nr <- exp(lnnr)	
	stockr$year <- y1:y2
	stockr$meanz <- apply(zmortr[,zbarage.1:zbarage.2], 1, mean)
	stockr$rec <- exp(rr[na:length(rr)])
	stockr$ssb <- apply(nr * swr * matr, 1, sum)	
	stockr$tsb <- apply(nr * swr, 1, sum)
	# Mean-standardise rec, ssb, tsb using full time-series means
	#mean.rec <- mean(stock$rec)
	#mean.ssb <- mean(stock$ssb)
	#mean.tsb <- mean(stock$tsb)
	#stockr$rec <- stockr$rec / mean.rec
	#stockr$ssb <- stockr$ssb / mean.ssb
	#stockr$tsb <- stockr$tsb / mean.tsb

	retro.stock[[kk]] <- stockr
}

# Plot retros (this also needs to be rewritten as a part of the 
# plot.surbar function)

if (!pdf.plots) windows(width = 10, height = 7)
par(mfrow = c(2,2), mar = c(5,4,2,1)+0.1)

# Mean Z
x.meanz <- do.call(rbind, lapply(x.psim.stock, function(wk){wk$meanz}))
x.meanz.quantile <- array(NA, dim = c(dim(x.meanz)[2],5))
x.meanz.mean <- rep(NA, dim(x.meanz)[2])
y2 <- max(as.numeric(rownames(sw)))
rownames(x.meanz.quantile) <- y1:y2
colnames(x.meanz.quantile) <- c("5%","25%","50%","75%","95%")
for (i in 1:dim(x.meanz)[2])
{
	x.meanz.quantile[i,] <- quantile(x.meanz[,i], c(0.05, 0.25, 0.5, 0.75, 0.95))
	x.meanz.mean[i] <- mean(x.meanz[,i])
}
plot(y1:y2, x.meanz.quantile[,3], type = "n", lty = 1,
	xlab = "Year", ylab = "Mean Z", 
	ylim = c(min(0, x.meanz.quantile), max(x.meanz.quantile)))
polygon(x = c(y1:y2, rev(y1:y2)), 
	y = c(x.meanz.quantile[,1], rev(x.meanz.quantile[,5])), density = -1, 
	col = "lightgrey", lty = 0)
lines(y1:y2, x.meanz.quantile[,3], lty = 1, col = "black", lwd = 2)
points(y2-1, x.meanz.quantile[ch(y2-1),3], pch = 16, col = "black")
for (k in 1:length(retro.stock))
{
	lines(retro.stock[[k]]$year, retro.stock[[k]]$meanz, lty = 1, col = "red")
	ly.less <- max(retro.stock[[k]]$year)-1 
	points(ly.less, retro.stock[[k]][retro.stock[[k]]$year == ly.less,]$meanz,
		pch = 16, col = "red")
}

# SSB
x.ssb <- do.call(rbind, lapply(x.psim.stock, function(wk){wk$ssb}))
x.ssb.quantile <- array(NA, dim = c(dim(x.ssb)[2],5))
x.ssb.mean <- rep(NA, dim(x.ssb)[2])
y2 <- max(as.numeric(rownames(sw)))
rownames(x.ssb.quantile) <- y1:y2
colnames(x.ssb.quantile) <- c("5%","25%","50%","75%","95%")
for (i in 1:dim(x.ssb)[2])
{
	x.ssb.quantile[i,] <- quantile(x.ssb[,i], c(0.05, 0.25, 0.5, 0.75, 0.95))
	x.ssb.mean[i] <- mean(x.ssb[,i])
}
plot(y1:y2, x.ssb.quantile[,3], type = "n", lty = 1,
	xlab = "Year", ylab = "SSB", 
	ylim = c(min(0, x.ssb.quantile), max(x.ssb.quantile)))
polygon(x = c(y1:y2, rev(y1:y2)), 
	y = c(x.ssb.quantile[,1], rev(x.ssb.quantile[,5])), density = -1, 
	col = "lightgrey", lty = 0)
lines(y1:y2, x.ssb.quantile[,3], lty = 1, col = "black", lwd = 2)
for (k in 1:length(retro.stock))
{
	lines(retro.stock[[k]]$year, retro.stock[[k]]$ssb, lty = 1, col = "red")
}

# TSB
x.tsb <- do.call(rbind, lapply(x.psim.stock, function(wk){wk$tsb}))
x.tsb.quantile <- array(NA, dim = c(dim(x.tsb)[2],5))
x.tsb.mean <- rep(NA, dim(x.tsb)[2])
y2 <- max(as.numeric(rownames(sw)))
rownames(x.tsb.quantile) <- y1:y2
colnames(x.tsb.quantile) <- c("5%","25%","50%","75%","95%")
for (i in 1:dim(x.tsb)[2])
{
	x.tsb.quantile[i,] <- quantile(x.tsb[,i], c(0.05, 0.25, 0.5, 0.75, 0.95))
	x.tsb.mean[i] <- mean(x.tsb[,i])
}
plot(y1:y2, x.tsb.quantile[,3], type = "n", lty = 1,
	xlab = "Year", ylab = "TSB", 
	ylim = c(min(0, x.tsb.quantile), max(x.tsb.quantile)))
polygon(x = c(y1:y2, rev(y1:y2)), 
	y = c(x.tsb.quantile[,1], rev(x.tsb.quantile[,5])), density = -1, 
	col = "lightgrey", lty = 0)
lines(y1:y2, x.tsb.quantile[,3], lty = 1, col = "black", lwd = 2)
for (k in 1:length(retro.stock))
{
	lines(retro.stock[[k]]$year, retro.stock[[k]]$tsb, lty = 1, col = "red")
}

# Recruitment
x.rec <- do.call(rbind, lapply(x.psim.stock, function(wk){wk$rec}))
x.rec.quantile <- array(NA, dim = c(dim(x.rec)[2],5))
x.rec.mean <- rep(NA, dim(x.rec)[2])
y2 <- max(as.numeric(rownames(sw)))
rownames(x.rec.quantile) <- y1:y2
colnames(x.rec.quantile) <- c("5%","25%","50%","75%","95%")
for (i in 1:dim(x.rec)[2])
{
	x.rec.quantile[i,] <- quantile(x.rec[,i], c(0.05, 0.25, 0.5, 0.75, 0.95))
	x.rec.mean[i] <- mean(x.rec[,i])
}
plot(y1:y2, x.rec.quantile[,3], type = "n", lty = 1,
	xlab = "Year", ylab = "Recruitment", 
	ylim = c(min(0, x.rec.quantile), max(x.rec.quantile)))
polygon(x = c(y1:y2, rev(y1:y2)), 
	y = c(x.rec.quantile[,1], rev(x.rec.quantile[,5])), density = -1, 
	col = "lightgrey", lty = 0)
lines(y1:y2, x.rec.quantile[,3], lty = 1, col = "black", lwd = 2)
for (k in 1:length(retro.stock))
{
	lines(retro.stock[[k]]$year, retro.stock[[k]]$rec, lty = 1, col = "red")
}
mtext("Retrospective analysis", outer = TRUE, side = 3, line = -1)

#####################################################################################

# Forecast - to be done!

#####################################################################################

# Stop timer and report
stop.time <- Sys.time()
cat(paste("\nRun complete in:", round(stop.time - start.time, 3), "s\n"))

# Close output
if (pdf.plots)
{
	dev.off()
}



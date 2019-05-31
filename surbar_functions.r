#############################################
# SURBA R implementation: functions         #
#                                           #
# Coby Needle, Marine Laboratory            #
#                                           #
# 14 December 2010                          #
# Modifications:                            #
# 31 May 2019 - modified for GitHub         #
# repository and TAF                        #
#############################################

# AIC function for nls.lm (Katherine Mullen)
logLik.nls.lm <- function(object, REML = FALSE, ...)
{
	res <- object$fvec
	N <- length(res)
	val <- -N * (log(2 * pi) + 1 - log(N) + log(sum(res^2)))/2
	## the formula here corresponds to estimating sigma^2.
	attr(val, "df") <- 1L + length(coef(object))
	attr(val, "nobs") <- attr(val, "nall") <- N
	class(val) <- "logLik"
	val
}
AIC.nls.lm <- function(object, REML = FALSE, ...)
{
	res <- object$fvec
	N <- length(res)
	val <- -N * (log(2 * pi) + 1 - log(N) + log(sum(res^2)))/2
	df <- 1L + length(coef(object))
	aic <- 2*N - 2*val 
	class(aic) <- "AIC"
	aic
}

#############################################
# FUNCTION: ch
#   
# Purpose: shorthand for as.character
# Created by: Coby Needle
# Date: 14/12/2010
# Modifications:
# Uses: 
#############################################

ch <- function(x){as.character(x)}

#############################################
# FUNCTION: read.vpa.file
#   
# Purpose: to read in Lowestoft VPA data
# Created by: Coby Needle
# Date: 14/12/2010
# Modifications:
# Uses: 
#############################################

read.vpa.file <- function(filename)
{
    wk.y <- scan(filename, skip = 2, nlines = 1, quiet = TRUE)
    wk.a <- scan(filename, skip = 3, nlines = 1, quiet = TRUE)
    wk.tab <- read.delim(filename, header = FALSE, sep = "", skip = 5)
    list(y = wk.y, a = wk.a, tab = wk.tab)
}

#############################################
# FUNCTION: read.survey.file
#   
# Purpose: to read in Lowestoft-format survey data
# Created by: Coby Needle
# Date: 14/12/2010
# Modifications:
# Uses: 
#############################################

read.survey.file <- function(filename)
{
    wk.n <- scan(filename, skip = 1, nlines = 1, quiet = TRUE) - 100
    wk.idx <- vector("list", length = wk.n)
    wk.start <- 3
    for (wk.k in 1:wk.n)
    {
            wk.idx[[wk.k]]$name <- paste(scan(filename, skip = wk.start - 1, nlines = 1, 
                what = character(0), quiet = TRUE), collapse = " ")
            wk.temp <- scan(filename, skip = wk.start, nlines = 1, quiet = TRUE)
            wk.idx[[wk.k]]$y1 <- wk.temp[1]
            wk.idx[[wk.k]]$y2 <- wk.temp[2]
            wk.idx[[wk.k]]$ny <- wk.temp[2] - wk.temp[1] + 1
            wk.temp <- scan(filename, skip = wk.start + 1, nlines = 1, quiet = TRUE)
            wk.idx[[wk.k]]$rho <- 0.5 * (wk.temp[4] + wk.temp[3])
            wk.temp <- scan(filename, skip = wk.start + 2, nlines = 1, quiet = TRUE)
            wk.idx[[wk.k]]$a1 <- wk.temp[1]
            wk.idx[[wk.k]]$a2 <- wk.temp[2]
            wk.idx[[wk.k]]$na <- wk.temp[2] - wk.temp[1] + 1
            wk.idx[[wk.k]]$tab <- read.table(filename, skip = wk.start + 3, nrows = wk.idx[[wk.k]]$ny)
            wk.temp <- wk.idx[[wk.k]]$tab[,2:(wk.idx[[wk.k]]$na + 1)] 
            wk.effort <- wk.idx[[wk.k]]$tab[,1] 
            wk.idx[[wk.k]]$tab <- data.frame(wk.temp / wk.effort)
            names(wk.idx[[wk.k]]$tab) <- wk.idx[[wk.k]]$a1:wk.idx[[wk.k]]$a2
            rownames(wk.idx[[wk.k]]$tab) <- wk.idx[[wk.k]]$y1:wk.idx[[wk.k]]$y2
            wk.start <- wk.start + 4 + wk.idx[[wk.k]]$ny
    }
    list(n = wk.n, idx = wk.idx)
}

#############################################
# FUNCTION: trim.data
#   
# Purpose: to trim survey indices to consistent year and age ranges
# Created by: Coby Needle
# Date: 14/12/2010
# Modifications:
# Uses: 
#############################################

trim.data <- function(wk.x, wk.idx, wk.qval, wk.wt, wk.sw, wk.mat)
{
    nk <- length(wk.x)
    a1 <- 100; a2 <- -100; y1 <- 9999; y2 <- -9999
    for (k in 1:nk)
    {
        ta1 <- wk.idx[[k]]$a1
        ta2 <- wk.idx[[k]]$a2
        ty1 <- wk.idx[[k]]$y1
        ty2 <- wk.idx[[k]]$y2
        if (ta1 < a1) a1 <- ta1 
        if (ta2 > a2) a2 <- ta2 
        if (ty1 < y1) y1 <- ty1 
        if (ty2 > y2) y2 <- ty2
    }
    a <- ch(a1:a2)
    y <- ch(y1:y2)
    for (k in 1:nk)
    {
        wk.x[[k]] <- wk.x[[k]][ch(y),ch(a)]
        wk.qval[[k]] <- wk.qval[[k]][ch(y),ch(a)]
        wk.wt[[k]] <- wk.wt[[k]][ch(y),ch(a)]
    }
    wk.sw <- wk.sw[ch(y),ch(a)]
    wk.mat <- wk.mat[ch(y),ch(a)]
    list(x = wk.x, qval = wk.qval, wt = wk.wt, sw = wk.sw, mat = wk.mat)
}

#############################################
# FUNCTION: trim.retro.data
#   
# Purpose: remove years from the end of each time series
# Created by: Coby Needle
# Date: 29/01/2011
# Modifications:
# Uses: 
#############################################

trim.retro.data <- function(wk.x, wk.qval, wk.wt, wk.sw, wk.mat)
{
    nk <- length(wk.x)
    y1 <- 9999; y2 <- -9999
    for (k in 1:nk)
    {
        ty1 <- min(row.names(wk.x[[k]]))
        ty2 <- max(row.names(wk.x[[k]]))
        if (ty1 < y1) y1 <- ty1 
        if (ty2 > y2) y2 <- ty2
    }
    y <- ch(y1:y2)
    for (k in 1:nk)
    {
        wk.x[[k]] <- wk.x[[k]][ch(y),]
        wk.qval[[k]] <- wk.qval[[k]][ch(y),]
        wk.wt[[k]] <- wk.wt[[k]][ch(y),]
    }
    wk.sw <- wk.sw[ch(y),]
    wk.mat <- wk.mat[ch(y),]
    list(x = wk.x, qval = wk.qval, wt = wk.wt, sw = wk.sw, mat = wk.mat)
}

#############################################
# FUNCTION: surbar
#   
# Purpose: SURBAR optimisation function
# Created by: Coby Needle
# Date: 14/12/2010
# Modifications:
# WKHAD 2014 - added an age-effect smoother
# Uses: 
#############################################

surbar <- function(wk.p)
{
    # Data x (list), qval, wt, rho.mod, lambda, spawn.time are assumed global
    
    wk.nk <- numk
    wk.na <- na
    wk.ny <- ny

    wk.s <- rep(NA, length = wk.na)
    wk.s[1:(ref.age-1)] <- wk.p[1:(ref.age-1)]
    wk.s[ref.age] <- 1.0
    wk.s[(ref.age+1):(wk.na-1)] <- wk.p[ref.age:(wk.na-2)]
    wk.s[wk.na] <- wk.s[wk.na-1]

    wk.f <- rep(NA, length = wk.ny)
    wk.f[1:(wk.ny-1)] <- wk.p[(wk.na-1):(wk.na+wk.ny-3)]
    wk.f[wk.ny] <- mean(wk.f[(wk.ny-3):(wk.ny-1)])

    wk.r <- rep(NA, length = wk.na + wk.ny - 1)
    wk.r[1:(wk.na+wk.ny-1)] <- wk.p[(wk.na+wk.ny-2):length(wk.p)]

    # Total mortality
    wk.z <- wk.f %o% wk.s

    # Abundance
    wk.n <- array(NA, dim = dim(wk.z))
    wk.n[1,] <- rev(wk.r[1:wk.na])
    wk.n[2:wk.ny,1] <- wk.r[(wk.na+1):length(wk.r)]

    vecs <- array(NA, dim = c(wk.na*wk.ny, 4))
    vecs[,1] <- matrix(data = wk.n, nrow = wk.na * wk.ny, ncol = 1)
    vecs[,2] <- rep(1:wk.na, each = wk.ny)
    vecs[,3] <- rep(y1:(y1 + wk.ny - 1), wk.na)
    vecs[,4] <- vecs[,3] - vecs[,2]

    cz.list <- tapply(wk.z, vecs[,4], cumsum)

    vecs.list <- lapply(levels(as.factor(vecs[,4])), function(wk){
        temp <- vecs[vecs[,4] == wk,]
        temp.rep <- dim(temp)[1]
        if (!is.null(temp.rep)) 
        {
            temp[,1] <- rep(temp[1], temp.rep)
        }
        temp
    })

    vecs.list <- lapply(vecs.list, function(wk){
        temp.rep <- dim(wk)[1]
        if (!is.null(temp.rep))
        {
            wk.zz <- unlist(cz.list[as.character(wk[1,4])])
            wk <- cbind(wk, wk.zz)
            wk.a <- dim(wk)[1]
            # lnN(a,y) <- lnN(a-1,y-1) - z(a-1,y-1)
            wk[2:wk.a,1] <- wk[1:(wk.a-1),1] - wk[1:(wk.a-1),5]
        }else
        {
            wk.zz <- unlist(cz.list[as.character(wk[4])])
            wk <- c(wk, wk.zz)
        }       
        wk
    })

    vecs.table <- do.call(rbind, vecs.list)
    vecs.table <- vecs.table[order(vecs.table[,3], vecs.table[,2]),]

    new.n <- matrix(vecs.table[,1], nrow = wk.ny, ncol = wk.na, byrow = TRUE)
    wk.n <- exp(new.n)

    # Fitted survey indices I.hat and
    # back-transformed observed survey indices I.dash.star
    # - now transformed to spawning time, not start of year - NEW CODE ################

    i.hat <- vector("list", length = wk.nk)
    i.dash.star <- vector("list", length = wk.nk)
    for (k in 1:wk.nk)
    {
        i.hat[[k]] <- wk.n * array(unlist(qval[[k]]), dim = dim(wk.n))
        i.dash.star[[k]] <- array(unlist(x[[k]] * exp(wk.z * (spawn.time - rho[k]))), dim = dim(wk.n))
    }

    # Survey log residuals
    res1 <- vector("list", length = wk.nk)
    out0 <- vector("list", length = wk.nk)
    for (k in 1:wk.nk)
    {
        res1[[k]] <- sqrt(wt[[k]]) * (log(i.dash.star[[k]]) - log(i.hat[[k]]))
        out0[[k]] <- array(NA, dim = c(wk.na * wk.ny, 1))
        out0[[k]][1:(wk.na*wk.ny),1] <- unlist(res1[[k]])
        out0[[k]] <- as.vector(na.exclude(out0[[k]]))
    }
    out1 <- as.vector(unlist(out0))

    # Lambda year-effect smoother
    f1 <- wk.f[1:(wk.ny-2)]
    f2 <- wk.f[2:(wk.ny-1)]
    res2 <- sqrt(lambda) * (f1 - f2)
    out2 <- as.vector(res2)

    # Lambda age-effect smoother
    s1 <- wk.s[1:(wk.na-2)]
    s2 <- wk.s[2:(wk.na-1)]
    res3 <- sqrt(lambda) * (s1 - s2)
    out3 <- as.vector(res3)

    # Overall SSQ
    as.numeric(c(out1,out2,out3))
}

#############################################
# FUNCTION: surbar.retro
#   
# Purpose: SURBAR retro optimisation function
# Created by: Coby Needle
# Date: 29/01/2011
# Modifications:
# Uses: 
#############################################

surbar.retro <- function(wk.p)
{
    # Data xr (list), qvalr, wtr, rho, lambda are assumed global
    
    wk.nk <- numk
    wk.na <- na
    wk.ny <- ny

    wk.s <- rep(NA, length = wk.na)
    wk.s[1:(ref.age-1)] <- wk.p[1:(ref.age-1)]
    wk.s[ref.age] <- 1.0
    wk.s[(ref.age+1):(wk.na-1)] <- wk.p[ref.age:(wk.na-2)]
    wk.s[wk.na] <- wk.s[wk.na-1]

    wk.f <- rep(NA, length = wk.ny)
    wk.f[1:(wk.ny-1)] <- wk.p[(wk.na-1):(wk.na+wk.ny-3)]
    wk.f[wk.ny] <- mean(wk.f[(wk.ny-3):(wk.ny-1)])

    wk.r <- rep(NA, length = wk.na + wk.ny - 1)
    wk.r[1:(wk.na+wk.ny-1)] <- wk.p[(wk.na+wk.ny-2):length(wk.p)]

    # Total mortality
    wk.z <- wk.f %o% wk.s

    # Abundance
    wk.n <- array(NA, dim = dim(wk.z))
    wk.n[1,] <- rev(wk.r[1:wk.na])
    wk.n[2:wk.ny,1] <- wk.r[(wk.na+1):length(wk.r)]

    vecs <- array(NA, dim = c(wk.na*wk.ny, 4))
    vecs[,1] <- matrix(data = wk.n, nrow = wk.na * wk.ny, ncol = 1)
    vecs[,2] <- rep(1:wk.na, each = wk.ny)
    vecs[,3] <- rep(y1:(y1 + wk.ny - 1), wk.na)
    vecs[,4] <- vecs[,3] - vecs[,2]

    cz.list <- tapply(wk.z, vecs[,4], cumsum)

    vecs.list <- lapply(levels(as.factor(vecs[,4])), function(wk){
        temp <- vecs[vecs[,4] == wk,]
        temp.rep <- dim(temp)[1]
        if (!is.null(temp.rep)) 
        {
            temp[,1] <- rep(temp[1], temp.rep)
        }
        temp
    })

    vecs.list <- lapply(vecs.list, function(wk){
        temp.rep <- dim(wk)[1]
        if (!is.null(temp.rep))
        {
            wk.zz <- unlist(cz.list[as.character(wk[1,4])])
            wk <- cbind(wk, wk.zz)
            wk.a <- dim(wk)[1]
            # lnN(a,y) <- lnN(a-1,y-1) - z(a-1,y-1)
            wk[2:wk.a,1] <- wk[1:(wk.a-1),1] - wk[1:(wk.a-1),5]
        }else
        {
            wk.zz <- unlist(cz.list[as.character(wk[4])])
            wk <- c(wk, wk.zz)
        }       
        wk
    })

    vecs.table <- do.call(rbind, vecs.list)
    vecs.table <- vecs.table[order(vecs.table[,3], vecs.table[,2]),]

    new.n <- matrix(vecs.table[,1], nrow = wk.ny, ncol = wk.na, byrow = TRUE)
    wk.n <- exp(new.n)

    # Fitted survey indices I.hat and
    # back-transformed observed survey indices I.dash.star
    # - now transformed to spawning time, not start of year - NEW CODE ################

    i.hat <- vector("list", length = wk.nk)
    i.dash.star <- vector("list", length = wk.nk)
    for (k in 1:wk.nk)
    {
        i.hat[[k]] <- wk.n * array(unlist(qvalr[[k]]), dim = dim(wk.n))
        i.dash.star[[k]] <- array(unlist(xr[[k]] * exp(wk.z * (spawn.time - rho[k]))), dim = dim(wk.n))
    }
    
    # Survey log residuals
    res1 <- vector("list", length = wk.nk)
    out0 <- vector("list", length = wk.nk)
    for (k in 1:wk.nk)
    {
        res1[[k]] <- sqrt(wtr[[k]]) * (log(i.dash.star[[k]]) - log(i.hat[[k]]))
        out0[[k]] <- array(NA, dim = c(wk.na * wk.ny, 1))
        out0[[k]][1:(wk.na*wk.ny),1] <- unlist(res1[[k]])
        out0[[k]] <- as.vector(na.exclude(out0[[k]]))
    }
    out1 <- as.vector(unlist(out0))

    # Lambda year-effect smoother
    f1 <- wk.f[1:(wk.ny-2)]
    f2 <- wk.f[2:(wk.ny-1)]
    res2 <- sqrt(lambda) * (f1 - f2)
    out2 <- as.vector(res2)

    # Lambda age-effect smoother
    s1 <- wk.s[1:(wk.na-2)]
    s2 <- wk.s[2:(wk.na-1)]
    res3 <- sqrt(lambda) * (s1 - s2)
    out3 <- as.vector(res3)

    # Overall SSQ
    as.numeric(c(out1,out2,out3))
}

#############################################
# FUNCTION: plot.index.corr
#   
# Purpose: bivariate scatterplot by age
# Created by: Colin Millar
# Date: unknown
# Modifications: updated by Coby Needle to work
# with non-FLR objects
# Uses: 
#############################################

# Functions to plot bivariate correlations

empty <- function(x, y, groups, subscripts, panel.number, packet.number) 
{
	# Empty function
}

diag.panel.cm <- function(draw, axis.line.col, ...)
{
	diag.panel.splom(draw = F, axis.line.col = axis.line.col, ...)
}

main.panel.cm <- function(x, y, groups, subscripts, 
	panel.number, packet.number, .tol = 0.05) 
{
	panel.xyplot(x, y, pch = 19, col = grey(0.35), cex = 0.2)

	# Fit a linear model to the data
	wk.lm1 <- lm(y ~ x)
	wk.x1 <- 0:20 / 20
	wk.fit <- suppressWarnings(predict.lm(wk.lm1, 
		newdata = data.frame(x = wk.x1), se.fit = T))
	wk.y1 <- wk.fit$fit
	wk.yu <- wk.y1 + 2*wk.fit$se
	wk.yl <- wk.y1 - 2*wk.fit$se

	wk.sig <- identical(anova(wk.lm1)$"Pr(>F)"[1] < .tol, T)

	if (wk.sig) 
	{
		wk.line.f <- list(lwd = 2, lty = 1, col = "#000000")
		wk.line.ci <- list(lwd = 1, lty = 1, col = "red4")
	} else 
	{
		wk.line.f <- list(lwd = 1, lty = 1, col = "#0000FF")
		wk.line.ci <- list(lwd = 1, lty = 2, col = "#0000FF")
	}

	panel.lines(wk.x1, wk.y1, lwd = wk.line.f$lwd, lty = wk.line.f$lty, col = wk.line.f$col)
	panel.lines(wk.x1, wk.yu, lwd = wk.line.ci$lwd, lty = wk.line.ci$lty, col = wk.line.ci$col)
	panel.lines(wk.x1, wk.yl, lwd = wk.line.ci$lwd, lty = wk.line.ci$lty, col = wk.line.ci$col)

	# Draw in axes
	grid.draw(linesGrob(x = c(-0.1,1,1), y = c(0,0,1.1),
		gp = gpar(lwd = 2, col = grey(0.5)),
		default.units = "npc"))
}

panel.pairs.cm <- function(z, subscripts, panel.subscripts) 
{
	panel.pairs(z,
		panel.subscripts = panel.subscripts,
		subscripts       = subscripts,
		lower.panel      = empty,
		diag.panel       = diag.panel.cm,
		upper.panel      = main.panel.cm,
		axis.line.col    = "#FFFFFF",
)

}

centre.log <- function(wk.mat) 
{
	wk.mat[wk.mat <= 0] <- NA
	wk.mat <- log(wk.mat)
	apply(wk.mat, 2, function(wk.x) (wk.x-min(wk.x, na.rm = TRUE))/diff(range(wk.x, na.rm = TRUE)))
}

plot.index.corr <- function(wk.object, wk.type) 
{
	# par(bty="n")
	# trellis.par.set(box.rectangle = list(col = "white"))

	for (wk.i in seq(length(wk.object))) 
	{
		# Select one tuning fleet
		if (wk.type == "FLR")
		{
			wk.tune.mat <- t(wk.object[[wk.i]]@catch.n@.Data[,,1,1,1,1])
			wk.main <- wk.object[[wk.i]]@name
		} else
		{
			wk.tune.mat <- wk.object[[wk.i]]
			wk.main <- names(wk.object)
		}
		# Make cohort matrix
		wk.n <- dim(wk.tune.mat)[2]
		wk.cohort.mat <- matrix(NA, ncol = wk.n, nrow = dim(wk.tune.mat)[1] + wk.n - 1)
		colnames(wk.cohort.mat) <- colnames(wk.tune.mat)
		for (wk.j in 1:wk.n) 
		{
			wk.cohort.mat[,wk.j] <- c(rep(NA,wk.n-wk.j), 
				wk.tune.mat[,wk.j], rep(NA,wk.j-1))
		}
		print(splom(~centre.log(wk.cohort.mat), superpanel = panel.pairs.cm, 
			xlab = wk.main, col = "white"))
	}
}


trim.surv <- function(wk.obj.index, wk.max.yr) 
{
	for (wk.i in 1:length(wk.obj.index)) 
	{
		wk.minyear <- wk.obj.index[[wk.i]]@range["minyear"]
		wk.maxyear <- wk.obj.index[[wk.i]]@range["maxyear"]
		if (wk.maxyear > wk.max.yr) wk.obj.index[[wk.i]] <- trim(wk.obj.index[[wk.i]], year = wk.minyear:wk.max.yr)
	}
	wk.obj.index
}

#############################################
# FUNCTION: plot.surbar
#   
# Purpose: general function to plot SURBAR outputs
# Created by: Coby Needle
# Date: 14/12/2010
# Modifications:
# Uses: 
#############################################

plot.surbar <- function(x, stock, psim, psim.s, psim.f, psim.r, res, y1, y2, 
    a1, a2, type, pdf)
{
    # Determine number of indices
    numk <- length(res)

    # Define window
    if (type == "sum.line" | type == "sum.boxplot")
    {
        if (!pdf)
        {
            windows(width = 10, height = 7)
        }
        par(mfrow = c(2,2), mar = c(5,4,1,1)+0.1)
    } else if (type == "res.line" | type == "res.smooth")
    {
        # if (!pdf) windows(width = 10, height = 7)
        #config <- switch(numk,
        #    c(1,1), c(2,1), c(2,2), c(2,2), c(2,3), c(2,3), c(3,3), c(3,3), c(3,3))
	  config <- c(1,1)
        par(mfrow = config, mar = c(5,4,1,1)+0.1)
    } else if (type == "catch.curve" | type == "log.by.cohort")
    {
        if (!pdf)
        {
            windows(width = 10, height = 7)
        }
        config <- switch(numk,
            c(1,1), c(2,1), c(2,2), c(2,2), c(2,3), c(2,3), c(3,3), c(3,3), c(3,3))
        par(mfrow = config, mar = c(5,4,3,1)+0.1)
    } else if (type == "params")
    {
        if (!pdf)
        {
            windows(width = 10, height = 7)
        }
        par(mfrow = c(2,3), mar = c(5,4,1,1)+0.1)
    } else if (type == "age.scatterplot")
    {
        # windows defined within function
    }

    # Generate plots
    if (type == "sum.line")
    {
        # Mean Z
        x.meanz <- do.call(rbind, lapply(psim, function(wk){wk$meanz}))
        x.meanz.quantile <- array(NA, dim = c(dim(x.meanz)[2],5))
        x.meanz.mean <- rep(NA, dim(x.meanz)[2])
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
            y = c(x.meanz.quantile[,1], rev(x.meanz.quantile[,5])), density = -1, col = "grey",
            lty = 0)
        lines(y1:y2, x.meanz.quantile[,3], lty = 1, col = "black", lwd = 2)
        points(stock$year, x.meanz.mean, pch = 3, col = "red")
        points(stock$year, stock$meanz, pch = 16, col = "green")
        legend(legend = c("NLS estimate", "Bootstrap mean", "Bootstrap median", "90% CI"), x = "bottomleft",
            pch = c(16,3,-1,-1), lty = c(-1,-1,1,1), lwd = c(-1,-1,2,5), bty = "n", col = c("green","red","black", "grey"))

        # SSB
        x.ssb <- do.call(rbind, lapply(psim, function(wk){wk$ssb}))
        x.ssb.quantile <- array(NA, dim = c(dim(x.ssb)[2],5))
        x.ssb.mean <- rep(NA, dim(x.ssb)[2])
        rownames(x.ssb.quantile) <- y1:y2
        colnames(x.ssb.quantile) <- c("5%","25%","50%","75%","95%")
        for (i in 1:dim(x.ssb)[2])
        {
            x.ssb.quantile[i,] <- quantile(x.ssb[,i], c(0.05, 0.25, 0.5, 0.75, 0.95))
            x.ssb.mean[i] <- mean(x.ssb[,i])
        }
        plot(y1:y2, x.ssb.quantile[,3], type = "n", lty = 1,
            xlab = "Year", ylab = "SSB", ylim = c(min(0, x.ssb.quantile), max(x.ssb.quantile)))
        polygon(x = c(y1:y2, rev(y1:y2)), 
            y = c(x.ssb.quantile[,1], rev(x.ssb.quantile[,5])), density = -1, col = "grey",
            lty = 0)
        lines(y1:y2, x.ssb.quantile[,3], lty = 1, col = "black", lwd = 2)
        points(stock$year, x.ssb.mean, pch = 3, col = "red")
        points(stock$year, stock$ssb, pch = 16, col = "green")

        # TSB
        x.tsb <- do.call(rbind, lapply(x.psim.stock, function(wk){wk$tsb}))
        x.tsb.quantile <- array(NA, dim = c(dim(x.tsb)[2],5))
        x.tsb.mean <- rep(NA, dim(x.tsb)[2])
        rownames(x.tsb.quantile) <- y1:y2
        colnames(x.tsb.quantile) <- c("5%","25%","50%","75%","95%")
        for (i in 1:dim(x.tsb)[2])
        {
            x.tsb.quantile[i,] <- quantile(x.tsb[,i], c(0.05, 0.25, 0.5, 0.75, 0.95))
            x.tsb.mean[i] <- mean(x.tsb[,i])
        }
        plot(y1:y2, x.tsb.quantile[,3], type = "n", lty = 1,
            xlab = "Year", ylab = "Total biomass", ylim = c(min(0, x.tsb.quantile), max(x.tsb.quantile)))
        polygon(x = c(y1:y2, rev(y1:y2)), 
            y = c(x.tsb.quantile[,1], rev(x.tsb.quantile[,5])), density = -1, col = "grey",
            lty = 0)
        lines(y1:y2, x.tsb.quantile[,3], lty = 1, col = "black", lwd = 2)
        points(stock$year, x.tsb.mean, pch = 3, col = "red")
        points(stock$year, stock$tsb, pch = 16, col = "green")

        # Recruitment
        x.rec <- do.call(rbind, lapply(x.psim.stock, function(wk){wk$rec}))
        x.rec.quantile <- array(NA, dim = c(dim(x.rec)[2],5))
        x.rec.mean <- rep(NA, dim(x.rec)[2])
        rownames(x.rec.quantile) <- y1:y2
        colnames(x.rec.quantile) <- c("5%","25%","50%","75%","95%")
        for (i in 1:dim(x.rec)[2])
        {
            x.rec.quantile[i,] <- quantile(x.rec[,i], c(0.05, 0.25, 0.5, 0.75, 0.95))
            x.rec.mean[i] <- mean(x.rec[,i])
        }
        plot(y1:y2, x.rec.quantile[,3], type = "n", lty = 1,
            xlab = "Year", ylab = "Recruitment", ylim = c(min(0, x.rec.quantile), max(x.rec.quantile)))
        polygon(x = c(y1:y2, rev(y1:y2)), 
            y = c(x.rec.quantile[,1], rev(x.rec.quantile[,5])), density = -1, col = "grey",
            lty = 0)
        lines(y1:y2, x.rec.quantile[,3], lty = 1, col = "black", lwd = 2)
        points(stock$year, x.rec.mean, pch = 3, col = "red")
        points(stock$year, stock$rec, pch = 16, col = "green")
    } else if (type == "sum.boxplot")
    {
        # Mean Z
        x.meanz <- do.call(rbind, lapply(psim, function(wk){wk$meanz}))
        bp.meanz <- data.frame(x.meanz)
        names(bp.meanz) <- y1:y2
        bp.meanz <- stack(bp.meanz)
        boxplot(values ~ ind, data = bp.meanz, xlab = "Year", ylab = "Mean Z")

        # SSB
        x.ssb <- do.call(rbind, lapply(psim, function(wk){wk$ssb}))
        bp.ssb <- data.frame(x.ssb)
        names(bp.ssb) <- y1:y2
        bp.ssb <- stack(bp.ssb)
        boxplot(values ~ ind, data = bp.ssb, xlab = "Year", ylab = "SSB")

        # TSB
        x.tsb <- do.call(rbind, lapply(psim, function(wk){wk$tsb}))
        bp.tsb <- data.frame(x.tsb)
        names(bp.tsb) <- y1:y2
        bp.tsb <- stack(bp.tsb)
        boxplot(values ~ ind, data = bp.tsb, xlab = "Year", ylab = "Total biomass")

        # Recruitment
        x.rec <- do.call(rbind, lapply(psim, function(wk){wk$rec}))
        bp.rec <- data.frame(x.rec)
        names(bp.rec) <- y1:y2
        bp.rec <- stack(bp.rec)
        boxplot(values ~ ind, data = bp.rec, xlab = "Year", ylab = "Recruitment")
    } else if (type == "res.line")
    {
        ymin <- min(as.numeric(unlist(res)))
        ymax <- max(as.numeric(unlist(res)))
        for (k in 1:numk)
        {
	        if (!pdf) windows(width = 10, height = 7)
            y1.a <- as.numeric(rownames(res[[k]])[1])
            y2.a <- as.numeric(rev(rownames(res[[k]]))[1])
    
            plot (y1:y2, rep(0, ny), xlab = "Year", ylab = "Log residual",
                ylim = c(ymin, ymax), type = "n")
            for (i in 1:dim(res[[k]])[2])
            {
                lines(y1.a:y2.a, res[[k]][,i], col = i)
            }
            abline(h = 0, col = 1)
            legend(x = "topleft", legend = idx[[k]]$name, bty = "n", cex = 0.75)    
            legend(x = "bottomleft", legend = paste("Age", names(res[[k]])), lty = 1, col = 1:dim(res[[k]])[2],
                bty = "n", cex = 0.75, ncol = 2)
        }
    } else if (type == "res.smooth")
    {
        if (!pdf) windows(width = 10, height = 7)
        ymin <- min(as.numeric(unlist(res)))
        ymax <- max(as.numeric(unlist(res)))
        for (k in 1:numk)
        {
            y1.a <- as.numeric(rownames(res[[k]])[1])
            y2.a <- as.numeric(rev(rownames(res[[k]]))[1])
    
            plot (y1:y2, rep(0, ny), xlab = "Year", ylab = "Log residual",
                ylim = c(ymin, ymax), type = "n")
            for (i in 1:dim(res[[k]])[2])
            {
                points(y1.a:y2.a, res[[k]][,i], pch = 3, col = i)
                data.loess <- data.frame(x = y1.a:y2.a, y = res[[k]][,i])
                res.loess <- loess(y ~ x, data = data.loess, span = 1.0)
                res.loess.pred.x <- seq(y1.a, y2.a, length = 100)
                res.loess.pred <- predict(res.loess, newdata = res.loess.pred.x, se = TRUE)
                lines(res.loess.pred.x, res.loess.pred$fit, col = i, lty = 1)
            }
            abline(h = 0, col = 1)
            legend(x = "topleft", legend = idx[[k]]$name, bty = "n", cex = 0.75)    
            legend(x = "bottomleft", legend = paste("Age", names(res[[k]])), lty = 1, col = 1:dim(res[[k]])[2],
                bty = "n", cex = 0.75, ncol = 2)
        }
    } else if (type == "params")
    {
        # Boxplots: s
        bp.s <- data.frame(psim.s)
        names(bp.s) <- a1:a2
        bp.s <- stack(bp.s)
        boxplot(values ~ ind, data = bp.s, xlab = "Age", ylab = "s")

        # Boxplots: f
        bp.f <- data.frame(x.psim.f)
        names(bp.f) <- y1:y2
        bp.f <- stack(bp.f)
        boxplot(values ~ ind, data = bp.f, xlab = "Year", ylab = "f")

        # Boxplots: r
        bp.r <- data.frame(x.psim.r)
        names(bp.r) <- (y1 - na + 1 - a1):(y2 - a1)
        bp.r <- stack(bp.r)
        boxplot(values ~ ind, data = bp.r, xlab = "Cohort", ylab = "r")

        # Lineplots: s
        x.s.quantile <- array(NA, dim = c(dim(psim.s)[2],5))
        x.s.mean <- rep(NA, dim(psim.s)[2])
        rownames(x.s.quantile) <- a1:a2
        colnames(x.s.quantile) <- c("5%","25%","50%","75%","95%")
        for (i in 1:dim(psim.s)[2])
        {
            x.s.quantile[i,] <- quantile(x.psim.s[,i], c(0.05, 0.25, 0.5, 0.75, 0.95))
            x.s.mean[i] <- mean(x.psim.s[,i])
        }
        plot(a1:a2, x.s.quantile[,3], type = "n", lty = 1,
            xlab = "Age", ylab = "s", ylim = c(min(0, x.s.quantile), max(x.s.quantile)))
        polygon(x = c(a1:a2, rev(a1:a2)), 
            y = c(x.s.quantile[,1], rev(x.s.quantile[,5])), density = -1, col = "grey",
            lty = 0)
        lines(a1:a2, x.s.quantile[,3], lty = 1, col = "black", lwd = 2)
        points(a1:a2, x.s.mean, pch = 3, col = "red")
        points(a1:a2, s, pch = 16, col = "green")
        legend(legend = c("NLS estimate", "Bootstrap mean", "Bootstrap median", "90% CI"), x = "bottomright",
            pch = c(16,3,-1,-1), lty = c(-1,-1,1,1), lwd = c(-1,-1,2,5), bty = "n", col = c("green","red","black", "grey"))

        # Lineplots: f
        x.f.quantile <- array(NA, dim = c(dim(psim.f)[2],5))
        x.f.mean <- rep(NA, dim(psim.f)[2])
        rownames(x.f.quantile) <- y1:y2
        colnames(x.f.quantile) <- c("5%","25%","50%","75%","95%")
        for (i in 1:dim(psim.f)[2])
        {
            x.f.quantile[i,] <- quantile(x.psim.f[,i], c(0.05, 0.25, 0.5, 0.75, 0.95))
            x.f.mean[i] <- mean(psim.f[,i])
        }
        plot(y1:y2, x.f.quantile[,3], type = "n", lty = 1,
            xlab = "Year", ylab = "f", ylim = c(min(0, x.f.quantile), max(x.f.quantile)))
        polygon(x = c(y1:y2, rev(y1:y2)), 
            y = c(x.f.quantile[,1], rev(x.f.quantile[,5])), density = -1, col = "grey",
            lty = 0)
        lines(y1:y2, x.f.quantile[,3], lty = 1, col = "black", lwd = 2)
        points(y1:y2, x.f.mean, pch = 3, col = "red")
        points(y1:y2, f, pch = 16, col = "green")

        # Lineplots: r
        r1 <- y1 - na + 1 - a1
        r2 <- y2 - a1
        x.r.quantile <- array(NA, dim = c(dim(psim.r)[2],5))
        x.r.mean <- rep(NA, dim(psim.r)[2])
        rownames(x.r.quantile) <- r1:r2
        colnames(x.r.quantile) <- c("5%","25%","50%","75%","95%")
        for (i in 1:dim(psim.r)[2])
        {
            x.r.quantile[i,] <- quantile(x.psim.r[,i], c(0.05, 0.25, 0.5, 0.75, 0.95))
            x.r.mean[i] <- mean(psim.r[,i])
        }
        plot(r1:r2, x.r.quantile[,3], type = "n", lty = 1,
            xlab = "Cohort", ylab = "r", ylim = c(min(0, x.r.quantile), max(x.r.quantile)))
        polygon(x = c(r1:r2, rev(r1:r2)), 
            y = c(x.r.quantile[,1], rev(x.r.quantile[,5])), density = -1, col = "grey",
            lty = 0)
        lines(r1:r2, x.r.quantile[,3], lty = 1, col = "black", lwd = 2)
        points(r1:r2, x.r.mean, pch = 3, col = "red")
        points(r1:r2, r, pch = 16, col = "green")
        abline(v = y1 - a1 - 0.5, lty = 8, col = "blue")
    } else if (type == "catch.curve")
    {
        for (k in 1:numk)
        {
            cc.years <- as.numeric(row.names(x[[k]]))
            cc.ages <- as.numeric(names(x[[k]]))

            cc <- stack(x[[k]])
            cc$ind <- as.numeric(as.character(cc$ind))
            cc <- data.frame(cbind(cc, cc.years))
            names(cc) <- c("index","age","year")
            cc <- data.frame(cbind(cc, cc$year - cc$age))
            names(cc) <- c("index","age","year","cohort")

            cc.cohorts <- unique(cc$cohort)[order(unique(cc$cohort))]
            plot(x = 0, y = 0, 
                xlim = c(min(cc.years, na.rm = TRUE), max(cc.years, na.rm = TRUE)),
                ylim = c(min(log(x[[k]]), na.rm = TRUE), max(log(x[[k]]), na.rm = TRUE) + 0.25),
                type = "n", xlab = "Year", ylab = "Log survey index",
                main = names(x)[k])

            jj <- 0
            for (j in cc.cohorts)
            {
                jj <- jj + 1
                cc0 <- cc[cc$cohort == j,]
                lines(cc0$year, log(cc0$index), type = "l", lty = 1, pch = -1,
                    col = rainbow(n = length(cc.cohorts))[jj])
                text(x = cc0$year[1], y = log(cc0$index)[1] + 0.25, cc0$cohort[1], cex = 0.75)
            }
        }
    } else if (type == "log.by.cohort")
    {

        # Set plot limits to cover all cohorts in all series
        lbc.xmin <- min(do.call(rbind, lapply(x, function(wk){as.numeric(rownames(wk))}))) -
            max(do.call(rbind, lapply(x, function(wk){as.numeric(names(wk))})))
        lbc.xmax <- max(do.call(rbind, lapply(x, function(wk){as.numeric(rownames(wk))})))

        for (k in 1:numk)
        {

            lbc <- x[[k]]

            # Remove all NA rows and columns (CLUMSY!)
            lbc.keeprows <- rep(TRUE, length = dim(lbc)[1])
            lbc.keepcols <- rep(TRUE, length = dim(lbc)[2])
            for (ii in 1:dim(lbc)[1]) # Rows
            {
                row.test <- !is.na(lbc[ii,])
                if (length(row.test[!row.test]) == length(row.test))
                {
                    lbc.keeprows[ii] <- FALSE
                }
            }
            for (jj in 1:dim(lbc)[2]) # Columns
            {
                col.test <- !is.na(lbc[,jj])
                if (length(col.test[!col.test]) == length(col.test))
                {
                    lbc.keepcols[jj] <- FALSE
                }
            }
            lbc <- lbc[lbc.keeprows,lbc.keepcols]

            # Mean standardise by age
            for (kk in 1:dim(lbc)[2])
            {
                lbc[,kk] <- lbc[,kk] / mean(lbc[,kk])
            }

            # Take logs
            lbc <- log(lbc)

            lbc.ages <- as.numeric(names(lbc))
            lbc.years <- as.numeric(rownames(lbc))

            # Generate stacked dataframe
            lbc <- stack(lbc)
            lbc$ind <- as.numeric(as.character(lbc$ind))
            lbc <- data.frame(cbind(lbc, lbc.years))
            names(lbc) <- c("index","age","year")
            lbc <- data.frame(cbind(lbc, lbc$year - lbc$age))
            names(lbc) <- c("index","age","year","cohort")

            lbc.cohorts <- unique(lbc$cohort)[order(unique(lbc$cohort))]
            plot(x = 0, y = 0, 
                xlim = c(lbc.xmin, lbc.xmax),
                ylim = c(min(lbc$index, na.rm = TRUE), max(lbc$index, na.rm = TRUE)),
                type = "n", xlab = "Cohort", ylab = "Mean-std log survey index",
                main = names(x)[k])

            jj <- 0
            for (j in lbc.ages)
            {
                jj <- jj + 1
                lbc0 <- lbc[lbc$age == j,]
                lines(lbc0$cohort, lbc0$index, type = "l", lty = 1, 
                    col = rainbow(n = length(lbc.ages))[jj])
            }
            jj <- 0
            for (j in lbc.ages)
            {
                jj <- jj + 1
                lbc0 <- lbc[lbc$age == j,]
                points(x = lbc0$cohort[1], y = lbc0$index[1], pch = 16, cex = 2.5, col = "white")
                points(x = lbc0$cohort[1], y = lbc0$index[1], pch = 1, cex = 2.5, col = rainbow(n = length(lbc.ages))[jj])
                text(x = lbc0$cohort[1], y = lbc0$index[1], lbc0$age[1], cex = 0.75)
            }
        }
    } else if (type == "age.scatterplot")
    {
        for (k in 1:numk)
        {
            plot.index.corr(x[k], pdf)
        }
    }
}

test.function <- function()
{
    print("This is test function...")
}

#############################################
# FUNCTION: surbar.gen
#   
# Purpose: Generates abundance and survey index values
# given SURBAR parameters (surbar in reverse)
# Created by: Coby Needle
# Date: 11/09/2013
# Modifications:
# Uses: 
#############################################

surbar.gen <- function(wk.p)
{
    # Data x.surba (list), qval, wt, rho.sim, lambda are assumed global
    
    wk.nk <- nk
    wk.na <- na
    wk.ny <- ny

    # wk.s <- rep(NA, length = wk.na)
    # wk.s[1:(ref.age-1)] <- wk.p[1:(ref.age-1)]
    # wk.s[ref.age] <- 1.0
    # wk.s[(ref.age+1):(wk.na-1)] <- wk.p[ref.age:(wk.na-2)]
    # wk.s[wk.na] <- wk.s[wk.na-1]

    # wk.f <- rep(NA, length = wk.ny)
    # wk.f[1:(wk.ny-1)] <- wk.p[(wk.na-1):(wk.na+wk.ny-3)]
    # wk.f[wk.ny] <- mean(wk.f[(wk.ny-3):(wk.ny-1)])

    # wk.r <- rep(NA, length = wk.na + wk.ny - 1)
    # wk.r[1:(wk.na+wk.ny-1)] <- wk.p[(wk.na+wk.ny-2):length(wk.p)]

	wk.s <- wk.p[1:wk.na]
	wk.f <- wk.p[(wk.na+1):(wk.na+wk.ny)]
	wk.r <- wk.p[(wk.na+wk.ny+1):length(wk.p)]

    # Total mortality
    wk.z <- wk.f %o% wk.s

    # Abundance
    wk.n <- array(NA, dim = dim(wk.z))
    wk.n[1,] <- rev(wk.r[1:wk.na])
    wk.n[2:wk.ny,1] <- wk.r[(wk.na+1):length(wk.r)]

    vecs <- array(NA, dim = c(wk.na*wk.ny, 4))
    vecs[,1] <- matrix(data = wk.n, nrow = wk.na * wk.ny, ncol = 1)
    vecs[,2] <- rep(1:wk.na, each = wk.ny)
    vecs[,3] <- rep(y1:(y1 + wk.ny - 1), wk.na)
    vecs[,4] <- vecs[,3] - vecs[,2]

    cz.list <- tapply(wk.z, vecs[,4], cumsum)

    vecs.list <- lapply(levels(as.factor(vecs[,4])), function(wk){
        temp <- vecs[vecs[,4] == wk,]
        temp.rep <- dim(temp)[1]
        if (!is.null(temp.rep)) 
        {
            temp[,1] <- rep(temp[1], temp.rep)
        }
        temp
    })

    vecs.list <- lapply(vecs.list, function(wk){
        temp.rep <- dim(wk)[1]
        if (!is.null(temp.rep))
        {
            wk.zz <- unlist(cz.list[as.character(wk[1,4])])
            wk <- cbind(wk, wk.zz)
            wk.a <- dim(wk)[1]
            # lnN(a,y) <- lnN(a-1,y-1) - z(a-1,y-1)
            wk[2:wk.a,1] <- wk[1:(wk.a-1),1] - wk[1:(wk.a-1),5]
        }else
        {
            wk.zz <- unlist(cz.list[as.character(wk[4])])
            wk <- c(wk, wk.zz)
        }       
        wk
    })

    vecs.table <- do.call(rbind, vecs.list)
    vecs.table <- vecs.table[order(vecs.table[,3], vecs.table[,2]),]

    new.n <- matrix(vecs.table[,1], nrow = wk.ny, ncol = wk.na, byrow = TRUE)
    wk.n <- exp(new.n)

    # Observed and forward-transformed survey indices
    i.hat <- vector("list", length = wk.nk)
    i.dash.star <- vector("list", length = wk.nk)
    for (k in 1:wk.nk)
    {
        i.hat[[k]] <- wk.n * array(unlist(qval[[k]]), dim = dim(wk.n))
        i.dash.star[[k]] <- i.hat[[k]] * exp(wk.z * -rho.sim[k])
    }

	# Return abundance and survey index
	list(wk.n = wk.n, i.hat = i.dash.star, wk.z = wk.z)
    
}




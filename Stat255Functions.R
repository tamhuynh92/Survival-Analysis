##
#####	Collection of functions written by D. Gillen that are used
#####	throughout Stat 255
#####
#####	Author: D. Gillen
#####	Course: Stat 255, Winter 2016
##
##
#####
#####	ifelse1() :	Helper function that allows for 
#####							returning a multivariate response
#####							based upon a univariate logical comparison
#####
##
ifelse1 <- function (test, yes, no){
    if (test) yes
    else no
		}
		
##		
#####	pKM() : Function to pull estimated survival probabilities
#####			at specified times for a fitted KM estimate
#####
#####			fit		:	an object of type survfit()
#####			q		:	quantiles for which probabilities are to be computed
#####			xscale 	:	an optional value to rescale time
##
pKM <- function( fit, q, xscale=1 ){
	fit$time <- fit$time / xscale
	if( is.null( fit$strata ) ) fit$strata <- length(fit$surv)
	nStrata <- length( fit$strata )
	rslt <- vector( "list", length=nStrata )
	names( rslt ) <- ifelse1( nStrata==1, "Survival Probabilities", names( fit$strata ) )
	for( i in 1:nStrata ){
		rslt[[i]] <- as.data.frame( cbind( q, matrix(NA, nrow=length(q), ncol=3) ) )
		names( rslt[[i]] ) <- c( "time", "km.est", "lower", "upper" ) 
		strata.fit <-  as.data.frame( cbind( fit$time, fit$surv, fit$lower, fit$upper )[ (c(0,cumsum(fit$strata))[i]+1):(c(0,cumsum(fit$strata))[i]+fit$strata[i]), ] )
		names( strata.fit ) <- c("time", "surv", "lower", "upper" )
		for(j in 1:length(q)){
			if( sum(strata.fit$time > q[j]) > 1 ) rslt[[i]][j,2:4] <- strata.fit[ max( which(strata.fit$time <= q[j] ) ), 2:4 ]
			}
		}
	return(rslt)
	}
	
	
##
#####	qKM() : Function to estimate quantiles of the survival distribution
#####			using a fitted KM estimate
#####
#####			fit		:	an object of type survfit()
#####			q		:	probabilities for which quantiles are to be computed
#####			xscale 	:	an optional value to rescale time
##
qKM <- function( fit, p, xscale=1 ){
	fit$time <- fit$time / xscale
	if( is.null( fit$strata ) ) fit$strata <- length(fit$surv)
	nStrata <- length( fit$strata )
	rslt <- vector( "list", length=nStrata )
	names( rslt ) <- ifelse1( nStrata==1, "Survival Quantiles", names( fit$strata ) )
	for( i in 1:nStrata ){
		rslt[[i]] <- as.data.frame( cbind( p, matrix(NA, nrow=length(p), ncol=3) ) )
		names( rslt[[i]] ) <- c( "perc.surv", "km.est", "lower", "upper" ) 
		strata.fit <-  as.data.frame( cbind( fit$time, fit$surv, fit$lower, fit$upper )[ (c(0,cumsum(fit$strata))[i]+1):(c(0,cumsum(fit$strata))[i]+fit$strata[i]), ] )
		names( strata.fit ) <- c("time", "surv", "lower", "upper" )
		for(j in 1:length(p)){
			if( sum(strata.fit$surv<=p[j]) ){
				rslt[[i]][j,2] <- strata.fit$time[ min( which(strata.fit$surv <=p[j] ) ) ]
				rslt[[i]][j,3] <- ifelse( sum(strata.fit$lower>p[j], na.rm=TRUE) >= 1, strata.fit$time[ max(which(strata.fit$lower>p[j]))+1 ], 0 )
				rslt[[i]][j,4] <- ifelse( sum(strata.fit$upper<p[j], na.rm=TRUE) >= 1, strata.fit$time[ min(which(strata.fit$upper<p[j])) ], Inf )
				}
			}
		}
	return(rslt)
	}	
	
##
#####	kmPlot() : Function to plot KM estimates that allows for labeling
#####				number of subjects at risk and number failures at multiple time points
#####
#####			fit		:	an object of type survfit()
#####			xscale 	:	an optional value to rescale time
#####			labelTimes : time points on x-axis where statistics should be computed
#####			groupLabels : vector of names used to define/label groups
##
kmPlot <- function( fit, xscale=365.25, xlab="Time from study start (yrs)", ylab="Survival", lty=1, 
					labelTimes=seq(0,max(fit$time),length=4)/xscale, groupLabels, xmarOffset=0, ymarOffset=0, cex=1 ){

	##	Compute number of strata and make sure labels match
	if( is.null( fit$strata ) ) fit$strata <- fit$n
	nStrata <- length(fit$strata)
	strata.var <- rep( 1:nStrata, times=fit$strata )
	if( missing(groupLabels) ) groupLabels <- paste( "Group", 1:nStrata )
		else if( length(groupLabels) != nStrata ) stop( "Group labels must be same length as number of strata" )
	
	##	Plot KM estimates 
	oldpar <- par( "mar", "mgp", "las" )
	par( mar=c(14-xmarOffset,14-ymarOffset,4,4)+.1, mgp=c(3,1.5,0), las=1 )
	plot( fit, ylab=ylab, xlab=xlab, xscale=xscale, lty=1:nStrata, mark.time=F, las=1, axes=F )
	abline( v=0 )
	axis( 1, at=seq(0, 100, by=.5 ) )
	axis( 2, at=seq(.2, 1, by=.2 ) )
	
	##	Compute and plot the number of subjects at risk and cum. number of events at each label time
	totalRisk <- totalEvents <- matrix( rep( NA, nStrata*length(labelTimes) ), nrow=nStrata )
	for( i in 1:nStrata ){
			text( -.05, -.3-.05*i, paste( groupLabels[i], "     ", sep="" ), adj=1, cex=cex*.7 )	
			mtext( side=1, at=-.2, line=i+4, paste( groupLabels[i], "     ", sep="" ), adj=1, cex=cex*.7, las=1 )	
			tempTime <- fit$time[ strata.var==i ]
			tempNrisk <- fit$n.risk[ strata.var==i ]		
			tempNevents <- fit$n.event[ strata.var==i ]
			for( j in 1:length(labelTimes) ){
				totalRisk[i,j] <- ifelse( j<=length(tempNrisk), tempNrisk[ sum( tempTime <= labelTimes[j]*xscale ) ], 0 )
				if( max(tempTime)<labelTimes[j]*xscale ) totalRisk[i,j] <- 0
				if( is.na( totalRisk[i,j] ) ) totalRisk[i,j] <- tempNrisk[1]
				totalEvents[i,j] <- sum( tempNevents[ tempTime <= (labelTimes[j]*xscale) ] ) 
				mtext( side=1, at=labelTimes[j], line=i+4, paste( totalRisk[i,j], " (", totalEvents[i,j], ")", sep="" ), cex=cex*.7 )
				if( i==nStrata & nStrata > 1 ) mtext( side=1, at=labelTimes[j], line=i+6, paste( sum(totalRisk[,j], na.rm=TRUE), 
										" (", sum(totalEvents[,j], na.rm=TRUE), ")", sep="" ), cex=cex*.7 )  
				}
			}
	if( nStrata > 1 ) mtext( side=1, at=-.2, line=i+6, paste( "Total   ", "     ", sep="" ), adj=1, cex=cex*.7 )
	invisible( par(oldpar) )
}

##
#####	survtrend() : Function to compute the Tarone trend test
#####
#####			formula	:	Surv() response and covariate for trend test
#####			data 	:	dataset containing response and predictor
#####			print.table : if TRUE, prints observed and expected 
#####						  failures under H_0 in each group
##
survtrend <- function( formula, data, print.table=TRUE ){
	lrfit <- survdiff( formula, data=data )
	df <- length( lrfit$n ) - 1
	score <- coxph( formula, data=data )$score
	if( print.table ){
		oetable <- cbind( lrfit$n, lrfit$obs, lrfit$exp )
		colnames(oetable) <- c("N", "Observed", "Expected" )
		print( oetable )
		}
	cat( "\nLogrank Test : Chi(", df, ") = ", lrfit$chisq, ", p-value = ", 1-pchisq(lrfit$chisq, df), sep="" )
	cat( "\nTarone Test Trend : Chi(1) = ", score, ", p-value = ", 1-pchisq(score, 1), sep="" )
	}

	
##
#####	Function to estimate linear contrasts of coefficients from a coxph() fit
##
linContr.coxph <- function( model, contr.names, contr.coef=rep(1,length(contr.names)) ){
	beta.hat <- model$coef 
	cov.beta <- vcov( model )

	contr.index <- match( contr.names, dimnames(cov.beta)[[1]] ) 
	beta.hat <- beta.hat[ contr.index ]
	cov.beta <- cov.beta[ contr.index,contr.index ]
	est <- contr.coef %*% beta.hat
	se.est <- sqrt( contr.coef %*% cov.beta %*% contr.coef )
	zStat <- est / se.est
	pVal <- 2*pnorm( abs(zStat), lower.tail=FALSE )
	ci95.lo <- exp( est - qnorm(.975)*se.est )
	ci95.hi <- exp( est + qnorm(.975)*se.est )
	est <- exp( est )

	cat( "\nTest of H_0: exp( " )
	for( i in 1:(length( contr.names )) ){
		if( i < length( contr.names ) ) cat( contr.coef[i], "*", contr.names[i], " + ", sep="" )
			else cat( contr.coef[i], "*", contr.names[i], " ) = 1 :\n\n", sep="" )	
		}

	rslt <- data.frame( est, se.est, zStat, pVal, ci95.lo, ci95.hi )
	colnames( rslt )[1] <- "exp( Est )"
	round( rslt, 3 )
}


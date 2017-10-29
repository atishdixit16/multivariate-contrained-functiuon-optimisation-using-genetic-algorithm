randomTrials <- function(count = 100, varCount = 2, range = matrix(c(-1,-1,1,1),2,2)) {
        trialStack <- NULL
        for (i in 1:count) {
		trial <- NULL
		for (j in 1:varCount)
			trial <- c( trial, runif(1,range[j,1],range[j,2]) )
		trialStack <- rbind(trialStack,trial)
	}
        trialStack
}

inputFunction <- function(values) {
	values[1]**2 + (values[2]-1)**2
}

constraints <- function(values, epsilon=0.0001) {
	max ( (abs(values[2] - values[1]**2) - epsilon) , 0 )
}

fitness <- function(values, func=inputFunction, constraints=TRUE, constrFunc=constraints, penalty=1) {
	if (constraints==FALSE)
		return(inputFunction(values))
	else 
		return( inputFunction(values) + penalty*constraints(values) )
}

getFitenssStack <- function(trialStacks, func= fitness, inputFunc=inputFunction, constraints=TRUE, constrFunc=constraints, penalty=1) {
        fitnessStack <- NULL
        for (i in 1:nrow(trialStacks))
                fitnessStack <- c( fitnessStack, func(trialStacks[i,], inputFunc,constraints=TRUE , constrFunc, penalty) )
        fitnessStack
}

naturalSelection <- function( trialStacks , type = "tournament", func = fitness, inputFunc=inputFunction, constraints=TRUE ,constrFunc=constraints, penalty=1 ) {
        if ( type=="tournament" ) {
                fitnessStack <- getFitenssStack (trialStacks, func, inputFunc,constraints=TRUE,constrFunc , penalty)
                len <- length(fitnessStack)
                winners <- NULL
                for (round in 1:len) {
                        pair <- sample(len,2)
                        if ( fitnessStack[ pair[1] ] <= fitnessStack[ pair[2] ] )
                                winners <- rbind(winners,trialStacks[pair[1],])
                        else
                                winners <- rbind(winners,trialStacks[pair[2],])
                }
                return (winners)
        }
}

crossOver <- function(trialStacks, crossProb = 0.8 ,type = "wholeArithmetic") {
        if (type=="wholeArithmetic") {
                len <- nrow(trialStacks)
                digits <- ncol(trialStacks)
		alpha <- runif(1)
                for (round in 1:len) {
                        pair <- sample(len,2)
                        if (runif(1) < crossProb) {
                                temp <- trialStacks[ pair[1], ]*alpha + trialStacks[ pair[2], ]*(1-alpha)
                                trialStacks[ pair[1], ]  <- trialStacks[ pair[2], ]*alpha + trialStacks[ pair[1], ]*(1-alpha)
                                trialStacks[ pair[2], ]  <- temp
                        }
                }
                return(trialStacks)
        }
}

mutation <- function(trialStacks, mutate_prob=(1/nrow(trialStacks)) , type = "realMutation" , timestep, totalTimestep, beta=1) {
        if (type=="realMutation") {
                len <- nrow(trialStacks)
                digits <- ncol(trialStacks)
		mins <- apply(trialStacks,2,min)
		maxs <- apply(trialStacks,2,max)
                for (i in 1:len)
                        for (j in 1:digits)
                                if (runif(1) < mutate_prob)
					if (runif(1) < 0.5)
                                        	trialStacks[i,j] <- trialStacks[i,j] + ( maxs[j] - trialStacks[i,j] )*(1-runif(1)^(1-timestep/totalTimestep))^beta
					else
                                        	trialStacks[i,j] <- trialStacks[i,j] - ( trialStacks[i,j] - mins[j] )*(1-runif(1)^(1-timestep/totalTimestep))^beta
                return(trialStacks)
        }
}

GeneOpt <- function( inputFunc=inputFunction, fitnessFunc =fitness, iterations=50, crossProb = 0.8, mutate_prob=0.001, count=100, varCount = 2, range = matrix(c(-1,-1,1,1),2,2), typeNS = "tournament", constraints=TRUE, constrFunc=constraints, penalty=1 , typeCross = "wholeArithmetic", typeMutate = "realMutation", beta=1 ) {

        trialStacks <- randomTrials(count, varCount, range)
        for ( i in 1:iterations) {
                trialStacks <- naturalSelection ( trialStacks , type = typeNS, func = fitnessFunc, inputFunc=inputFunc, constraints=constraints , constrFunc=constrFunc, penalty )
                trialStacks <- crossOver( trialStacks, crossProb = crossProb ,type = typeCross )
                trialStacks <- mutation ( trialStacks, mutate_prob=mutate_prob , type = typeMutate , i, iterations, beta )
        }
	sols <- getFitenssStack (trialStacks, fitnessFunc, inputFunc,constraints=FALSE,constrFunc , penalty)
        min_f <- min( sols )
	arg_min_f = unique ( trialStacks[ which ( sols == min_f ),] )
	#output
	list( min = min_f , arg_min = arg_min_f)
}

Answer <- GeneOpt( inputFunc=inputFunction, fitnessFunc =fitness, iterations=1000, crossProb = 0.8, mutate_prob=0.001, count=100, varCount = 2, range = matrix(c(-1,-1,1,1),2,2), typeNS = "tournament", constraints=TRUE, constrFunc=constraints, penalty=1 , typeCross = "wholeArithmetic", typeMutate = "realMutation", beta=1 )

print(Answer)

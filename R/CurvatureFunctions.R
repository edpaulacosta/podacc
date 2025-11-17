#' @export
invertFunction <- function(x,y){
  
  newX = max(x)-x
  newY = max(y)-y
  
  return(cbind(newX,newY))
}

#' @export
normalize <- function(a){
  return ((a - min(a)) / (max(a) - min(a)))
}

#' @export
invertNormalizeValue <- function(value,a){
  
  newValue = (max(a) - min(a))*value+min(a)
  
  return(newValue)
}


#' @export
maxDiffStraightLine = function(x,y){
  
  xN = normalize(x)
  yN = normalize(y)
  
  yDifference = yN-xN
  
  return(x[which.max(abs(yDifference))])
  
}


#' @export
returnArrayKnees <- function(x,y,S=1,curve = 'concave',direction = 'increasing'){
  
  
  knee = NA
  
  #library(reticulate)
  
  originalX = x
  originalY = y
  
  #sci = reticulate::import("scipy.interpolate")
  #uspline = sci$interp1d(x, y)
  #Ds_y = uspline(x)
  
  uspline = smooth.spline(x,y)
  Ds_y = predict(uspline,x=x)
  
  
  
  x_normalized = normalize(x)
  y_normalized = normalize(Ds_y$y)
  
  mDiff = mean(x_normalized-y_normalized)
  
  #plot(x_normalized,y_normalized)
  
  # the difference curve (add 1 so that it is never negative)
  # Added this line
  if(length(mDiff)==0){
    return(NA)
  }
  
  if(mDiff>0){
    y_difference = x_normalized-y_normalized 
    
  }else{
    y_difference = y_normalized-x_normalized
  }
  
  x_difference = x_normalized
  
  
  #np = reticulate::import("numpy")
  #ss = reticulate::import("scipy.signal")
  #YDiff = np$array(y_difference)
  #maxima_indices = ss$argrelextrema(YDiff, np$greater)[[1]]
  
  
  maxima_indices = which(diff(sign(diff(y_difference)))==-2)
  maxima_indices = maxima_indices+1 #account for difference between python and R
  
  
  x_difference_maxima = x_difference[maxima_indices]
  y_difference_maxima = y_difference[maxima_indices]
  
  
  # local minima
  #minima_indices = ss$argrelextrema(YDiff, np$less)[[1]]
  minima_indices = which(diff(sign(diff(y_difference)))==2)
  minima_indices = minima_indices+1
  x_difference_minima = x_difference[minima_indices]
  y_difference_minima = y_difference[minima_indices]
  
  min=950
  max=1000
  min=1
  max=length(x_difference)
  y_differenceNew = y_difference[order(x_difference)]
  x_differenceNew = x_difference[order(x_difference)]
  
  
  newX = seq(min(x_normalized),max(x_normalized),by=0.1)
  newY = rep(y_difference_maxima,length(newX))
  
  
  # Sensitivity parameter S
  # smaller values detect knees quicker
  
  
  
  #Tmx = abs(y_difference_maxima) - (S * np$abs(np$mean(np$diff(x_normalized))))
  Tmx = abs(y_difference_maxima) - (S * abs(mean(diff(x_normalized))))
  
  # artificially place a local max at the last item in the x_difference array
  #maxima_indices = np$append(maxima_indices, length(x_difference))
  maxima_indices = c(maxima_indices, length(x_difference))
  
  #minima_indices = np$append(minima_indices, length(x_difference))
  minima_indices = c(minima_indices, length(x_difference))
  
  # placeholder for which threshold region i is located in.
  maxima_threshold_index = 1
  minima_threshold_index = 1
  
  
  all_knees = c("")
  all_norm_knees = c("")
  
  
  curve = 'convex'
  direction = 'decreasing'
  
  
  tmp=cbind(x,as.numeric(as.character(y_difference)))
  tmp=cbind(tmp,as.numeric(as.character(x_normalized)))
  
  
  thresholIdx = rep(0,length(y_difference))
  knee = rep("No",length(y_difference))
  
  threshold_indices = rep(0,length(y_difference))
  tmp=cbind(tmp,thresholIdx)
  tmp=cbind(tmp,knee)
  tmp=cbind(tmp,threshold_indices)
  
  colnames(tmp) = c("x","y_difference","x_normalized","thresholIdx","knee","threshold_indices")
  
  
  # traverse the difference curve
  for (idx in 1:length(x_difference)){
    #for (idx in 1:78){
    
    
    #if(x_difference[idx] == 1.0){
    #  break
    #}
    
    # values in difference curve are at or after a local maximum
    if (idx >= maxima_indices[maxima_threshold_index]){
      if(idx<length(x_difference)){	#allow the last point to help in the call 
        threshold = Tmx[maxima_threshold_index]
        threshold_index = idx
        maxima_threshold_index = maxima_threshold_index + 1
      }	
    }
    
    # values in difference curve are at or after a local minimum
    if (idx >= minima_indices[minima_threshold_index]){
      if(idx<length(x_difference)){	#allow the last point to help in the call        
        
        threshold = 0.0
        minima_threshold_index = minima_threshold_index + 1
      }
    }
    
    
    # Do not evaluate values in the difference curve before the first local maximum.
    if (idx < maxima_indices[1]){
      next
    }
    
    
    tmp[idx,"thresholIdx"]=threshold
    
    
    tmp[idx,"threshold_indices"]=threshold_index 
    
    if (abs(y_difference[idx]) < abs(threshold)){
      tmp[idx,"knee"]=x[threshold_index]
      knee = originalX[threshold_index]
      all_knees = c(all_knees,knee)
      norm_knee = x_normalized[threshold_index]
      all_norm_knees = c(all_norm_knees,norm_knee)
      
    }
    
  }
  
  all_knees = union(all_knees,all_knees)
  all_norm_knees = union(all_norm_knees,all_norm_knees)
  
  all_knees = setdiff(all_knees,"")
  all_norm_knees = setdiff(all_norm_knees,"")
  
  
  
  
  
  return(all_knees)
  
  
  
}

#' @export
findElbow <- function(bmdValues,accumulated,Log=TRUE,S=1){
  
  # Gives the option to look at log10 transformed curve or the non-transformed data
  if(Log){
    bmdValuesLOG = log(bmdValues,10)
  }else{
    bmdValuesLOG =bmdValues
  }	
  
  output=invertFunction(bmdValuesLOG,accumulated)
  bmdValuesInv=output[,1]
  accumulatedInv=output[,2]
  #plot(bmdValues,accumulated,log="x")
  #plot(bmdValuesInv,accumulatedInv,col="blue")
  
  
  newBmdValuesInv = bmdValuesInv
  newAccumulatedInv= accumulatedInv
  
  
  x=newBmdValuesInv[order(newBmdValuesInv)]
  y=newAccumulatedInv[order(newBmdValuesInv)]
  #plot(x,y,xlab="BMD", ylab="Accumulation")
  
  
  allKnees=returnArrayKnees(x,y)
  
  
  #plot(bmdValuesInv,accumulatedInv,col="blue")
  #abline(v=allKnees[length(allKnees)],col="green")
  
  
  
  if(length(allKnees)==0){
    return(list(NA,NA))
  }else{
    
    # Mapping back to the original function
    # Loop to print all elbows if it is necessary
    elbows = rep(NA,length(allKnees))
    accumulated_elbow = rep(NA,length(allKnees))
    for(k in 1:length(allKnees)){
      mask = (allKnees[k]==bmdValuesInv)
      if(Log){
        elbows[k] = 10^(bmdValuesLOG[mask])
        accumulated_elbow[k] = accumulated[mask]
      }else{
        elbows[k] = bmdValuesLOG[mask]
      }
      #abline(v=elbow,col="red")
    }
    #mask = (allKnees[length(allKnees)]==bmdValuesInv)
    #if(Log){
    #	elbow = 10^(bmdValuesLOG[mask])
    #}else{
    #	elbow = bmdValuesLOG[mask]
    #}
    
    output = list(elbows,accumulated_elbow)
    
  }
  
  return(output)
  
}









#' @export
findIntersectionLine <- function(x,y){
  
  #sci = reticulate::import("scipy.interpolate")
  #uspline = sci$interp1d(x, y)
  #Ds_y = uspline(x)
  
  uspline = smooth.spline(x,y)
  Ds_y = predict(uspline,x=x)
  
  
  
  x_normalized = normalize(x)
  y_normalized = normalize(Ds_y$y)
  
  #plot(x_normalized,y_normalized)
  #points(x_normalized,x_normalized,col="red")
  #points(x_normalized,y_normalized,col="yellow")
  
  diffLine = x_normalized-y_normalized
  diffLine=round(diffLine,2)
  
  diffLine = diffLine[!is.na(diffLine)]
  #plot(diffLine)
  
  
  if(length(diffLine)==0){
    return(-1)
  }
  
  # If curve is close to a line
  equal_zero = length(diffLine[diffLine<=0.01&diffLine>=0.01])
  if((equal_zero/length(diffLine))>0.95){  
    return(-1)
  }
  
  # Remove zeros and negative values at the beginning
  # count number of removed values
  flag=TRUE
  count = 0  
  if(length(diffLine)>1){
    if(diffLine[1]>0){
      flag=FALSE
    }else{
      diffLine=diffLine[-1]
      if(diffLine[1]<0){
        count = count+1
      }
    }
  }
  
  while(flag){		
    if(length(diffLine)>1){
      diffLine=diffLine[-1]
      if(diffLine[1]<0){
        count = count+1
      }
      if(diffLine[1]>0){
        flag=FALSE
      }
    }else{
      flag=FALSE
    }
    
  }
  
  if(count>20){ # the identity line is starting below the curve
    # we will try the next iterations to find a valid intersection point
    return(1)
  }
  
  # if there are only negtaive number we will not output any value
  if(length(diffLine)==1){
    posPointFinal = -1
  }else{  
    
    point = diffLine[diffLine<0]
    
    # if point == 0, it is the end of the curve.   
    if(length(point)==0){
      posPointFinal = length(x)
    }else{
      pointFinal = point[1]
      diffLine[diffLine<0]=rep(-1,length(point))	
      posPointFinal = which.min(diffLine)
      #abline(v=x_normalized[posPointFinal],col="blue")
    }
  }
  
  return(posPointFinal)
  
  
  
}


#' @export
chooseSmoothing = function(bmdValuesLOG,accumulated){ 
  
  #library(reticulate)
  #np = reticulate::import("numpy")
  
  pp=smooth.spline(as.numeric(bmdValuesLOG),as.numeric(accumulated),all.knots=TRUE,df=20)
  
  #intervalX = np$abs(np$mean(np$diff(bmdValuesLOG)))
  intervalX = abs(mean(diff(bmdValuesLOG)))
  
  
  newBmdValuesLOG = seq(from = min(bmdValuesLOG), to = max(bmdValuesLOG), by = intervalX/10)
  newY = predict(pp,x=newBmdValuesLOG)
  
  if(is.list(newY)){
    newY = newY$y
  }else{
    newY = as.numeric(predict(fit,x=newBmdValuesLOG))
  }
  
  max = newY[1]
  flag="pp"
  h=2
  flagLoop=TRUE	
  while(flagLoop){
    
    cur = newY[h]
    if(cur>max){
      max=cur
    }else{
      if((max-cur)>0.01){
        flag="scam"
        flagLoop=FALSE
      }
    }
    
    if(h==length(newY)){
      flagLoop=FALSE		
    }
    h=h+1
    
  }
  
  return(flag)
}


# Function to check if values decrease at least once
check_decrease <- function(arr) {
  any(diff(arr) < 0)
}

# Function to check for decrease and get previous value
get_previous_value <- function(arr) {
  diffs <- diff(arr)
  decrease_indices <- which(diffs < 0)
  
  for(i in 1:length(decrease_indices)){
    
    arr[decrease_indices[i]+1]=arr[decrease_indices[i]]
    
  }
  return(arr)
}



ensure_curve_only_increases = function(pp){
  
  while(check_decrease(pp$y)){
    pp$y=get_previous_value(pp$y)
  }
  
  return(pp) 
}







#' @export
performSmoothing = function(bmdValuesLOG,accumulated,dfreedom,choiceSmoothing="",scam.m=2){
  
  library("scam")
  
  #np = reticulate::import("numpy")
  
  if(length(bmdValuesLOG)<10){
    choiceSmoothing = "NO"
  }
  
  if(choiceSmoothing==""){
    choiceSmoothing = chooseSmoothing(bmdValuesLOG,accumulated)
  }
  
  intervalX = abs(mean(diff(bmdValuesLOG)))
  if(intervalX>0.02374724){
    intervalX=0.02374724
  }
  
  
  newBmdValuesLOG = seq(from = min(bmdValuesLOG), to = max(bmdValuesLOG), by = intervalX)
  
  if(choiceSmoothing=="NO"){
    newBmdValuesLOG = bmdValuesLOG
  }
  
  if(choiceSmoothing=="pp"){
    pp=smooth.spline(as.numeric(bmdValuesLOG),as.numeric(accumulated),all.knots=TRUE,df=dfreedom)	
    pp=ensure_curve_only_increases(pp)
    newY = predict(pp,x=newBmdValuesLOG)
    #plot(pp)
    #plot(newY)
    #plot(as.numeric(newBmdValuesLOG[,1]),as.numeric(newY$y[,1]))
    fit = pp
    accumulated = newY$y
    
  }
  if(choiceSmoothing=="scam"){
    
    #fit=try(scam(accumulated~s(bmdValuesLOG,k=20,bs="mpi",m=2),family=gaussian))
    
    m=scam.m
    
    fit=try(scam(accumulated~s(bmdValuesLOG,k=dfreedom,bs="mpi",m=m),family=gaussian(link="identity")))
    if(length(fit)==1){
      fit=scam(accumulated~s(bmdValuesLOG,bs="mpi",m=m),family=gaussian(link="identity"))
    }
    
    
    #plot(fit)
    bmdValuesLOG = as.data.frame(bmdValuesLOG)
    newBmdValuesLOG = as.data.frame(newBmdValuesLOG)
    names(newBmdValuesLOG) = names(bmdValuesLOG)
    newY = as.numeric(predict(fit,newdata=newBmdValuesLOG))
    #plot(newY)
    #plot(newBmdValuesLOG,newY)
    accumulated = newY
    newBmdValuesLOG = newBmdValuesLOG[,1]
  }
  if(choiceSmoothing=="new"){
    
    m=scam.m
    pp=smooth.spline(as.numeric(bmdValuesLOG),as.numeric(accumulated),all.knots=TRUE,df=dfreedom)
    pp=ensure_curve_only_increases(pp)
    #plot(pp)
    newY = predict(pp,x=newBmdValuesLOG)
    fit=try(scam(newY$y~s(newBmdValuesLOG,k=dfreedom,bs="mpi",m=m),family=gaussian(link="identity")))
    if(length(fit)==1){
      fit=scam(newY$y~s(newBmdValuesLOG,bs="mpi",m=m),family=gaussian(link="identity"))
    }
    newBmdValuesLOG = as.data.frame(newBmdValuesLOG)
    newY =try(predict(fit,newdata=newBmdValuesLOG))
    if(length(newY)==1){
      pp=smooth.spline(as.numeric(bmdValuesLOG),as.numeric(accumulated),all.knots=TRUE,df=dfreedom)	
      newY = predict(pp,x=newBmdValuesLOG)
      fit = pp
      accumulated = newY$y
      accumulated = accumulated[,1]
      newBmdValuesLOG = newBmdValuesLOG[,1]
    }else{
      #plot(newY)
      accumulated = as.numeric(newY)
      newBmdValuesLOG = newBmdValuesLOG[,1]
    }
  }
  
  if(choiceSmoothing!="NO"){
    pp=smooth.spline(as.numeric(newBmdValuesLOG),as.numeric(accumulated),all.knots=TRUE,df=dfreedom)
  }else{
    pp=smooth.spline(as.numeric(newBmdValuesLOG),as.numeric(accumulated),all.knots=TRUE)
  }
  pp=ensure_curve_only_increases(pp)
  
  mask = pp$y<1
  pp$y[mask]=1
  
  
  return(pp)
  
}


#' @export
extrapolateDots = function(bmdValuesLOG,accumulated,min_interval=5){
  # What does it do?
  # Inserts artifically dots to help in the smoothing process
  
  intervals_x = diff(bmdValuesLOG)
  
  if(max(intervals_x)>min_interval/100){
    #if(max(intervals_x)/median(intervals_x)>min_interval){
    
    #cat(max(intervals_x)/median(intervals_x))
    #cat("\n")
    
    index_max = which.max(intervals_x)
    
    new_bmd = (bmdValuesLOG[index_max+1]+bmdValuesLOG[index_max])/2
    #new_accumulated = (accumulated[index_max+1]+accumulated[index_max])/2
    new_accumulated = (accumulated[index_max])
    
    bmdValuesLOG = c(bmdValuesLOG,new_bmd)
    accumulated = c(accumulated,new_accumulated)
    
    accumulated = accumulated[order(bmdValuesLOG)]
    bmdValuesLOG = sort(bmdValuesLOG)
    
    #cat(length(accumulated))
    #cat("\n")
    
    output = extrapolateDots(bmdValuesLOG,accumulated,min_interval=min_interval)
    bmdValuesLOG = output[,1]
    accumulated = output[,2]
    
  }
  
  #plot(bmdValuesLOG,accumulated)
  
  return(cbind(bmdValuesLOG,accumulated))
  
}

#' @export
extrapolateDotsUpToLength = function(bmdValuesLOG,accumulated,length=40){
  
  if(length(bmdValuesLOG)<length){
    
    intervals_x = diff(bmdValuesLOG)
    
    index_max = which.max(intervals_x)
    
    new_bmd = (bmdValuesLOG[index_max+1]+bmdValuesLOG[index_max])/2
    new_accumulated = (accumulated[index_max+1]+accumulated[index_max])/2
    
    bmdValuesLOG = c(bmdValuesLOG,new_bmd)
    accumulated = c(accumulated,new_accumulated)
    
    accumulated = accumulated[order(bmdValuesLOG)]
    bmdValuesLOG = sort(bmdValuesLOG)
    
    output = extrapolateDotsUpToLength(bmdValuesLOG,accumulated,length)
    bmdValuesLOG = output[,1]
    accumulated = output[,2]
    
    
  }
  
  #plot(bmdValuesLOG,accumulated)
  
  return(cbind(bmdValuesLOG,accumulated))
  
}


#' @export
getSmoothedCurve = function(bmdValues,accumulated,limit_curve_up,startPos=1,choiceSmoothing="pp",dfreedom=20,scam.m=2){
  
  bmdValues = bmdValues[startPos:length(bmdValues)]
  accumulated = accumulated[startPos:length(accumulated)]
  
  accumulated_original = accumulated
  
  if(length(bmdValues[bmdValues<=limit_curve_up])<5){
    return(NA)
  }else{
    
    #plot(bmdValues,accumulated,log="x")
    
    bmdValuesLOG = as.numeric(log(bmdValues,10))
    
    # extrapolate dots to allow smoothing for extreme datasets
    # if we dont put limits the recursion will cause error
    if(length(bmdValuesLOG)<500){
      min_interval= 5
    }else{
      min_interval = 20
    }
    
    output = extrapolateDots(bmdValuesLOG,accumulated,min_interval = min_interval)
    bmdValuesLOG = output[,1]
    accumulated = output[,2]
    #plot(bmdValuesLOG,accumulated)
    #bmdValuesLOG2 = output[,1]
    #accumulated2 = output[,2]
    #plot(bmdValuesLOG2,accumulated2)
    
    
    output = extrapolateDotsUpToLength(bmdValuesLOG,accumulated,length=40)
    bmdValuesLOG = output[,1]
    accumulated = output[,2]
    
    
    if(dfreedom>length(bmdValuesLOG)){
      dfreedom=round(2*(length(bmdValuesLOG)/3))
      
    }else{
      dfreedom=min(max(round(length(bmdValuesLOG)/10),20),40)
    }
    
    #length(bmdValuesLOG)
    
    pp = performSmoothing(as.numeric(bmdValuesLOG),accumulated,dfreedom,choiceSmoothing=choiceSmoothing,scam.m=scam.m)
    
    #pp = performSmoothing(as.numeric(bmdValuesLOG),accumulated,dfreedom,choiceSmoothing="pp",scam.m=scam.m)
    
    #plot(pp)
    #points(pp$x,pp$y,col="yellow")	
    #points(log(bmdValues,10),accumulated_original,col="red")
    #plot(log(bmdValues,10),accumulated,col="red")
    #plot(bmdValues,accumulated,log="x")
    
    mask = pp$x<=log(limit_curve_up,10)
    pp$x = pp$x[mask]
    pp$y = pp$y[mask]
    
    mask=pp$y<pp$y[1]
    pp$y[mask] = pp$y[1]
    
    ppOriginal = pp
    
    # Finding intersection between the initial point and the limited portion of the curve
    if(length(pp$y)>0){
      posPointFinal = findIntersectionLine(pp$x,pp$y)
      firstIntersection = posPointFinal
      # If it is the end of the curve (second IF condition. Then no need to further test)
      if((posPointFinal!=-1)&(posPointFinal!=length(pp$x))){
        # Limiting the size of the curve of 5 to look for intersection
        while(firstIntersection<length(pp$x)&length(pp$x)>5){
          pp$x = pp$x[1:(length(pp$x)-1)]
          pp$y = pp$y[1:(length(pp$y)-1)]
          posPointFinal = findIntersectionLine(pp$x,pp$y)
          if(posPointFinal>firstIntersection){
            firstIntersection = posPointFinal
          }
          #abline(v=pp$x[length(pp$x)],col="red")
          #abline(v=pp$x[posPointFinal],col="red")
          #cat(posPointFinal)
          #cat("\t")
          #cat(length(pp$x))
          #cat("\n")
        }		
        
        # If values ends up being 1, it means no intersection was ever found because the unity line was under the curve
        if(firstIntersection==1){
          posPointFinal=-1
        }else{
          
          posPointFinal = firstIntersection
          
          pp$x = ppOriginal$x[1:posPointFinal]
          pp$y = ppOriginal$y[1:posPointFinal]
          
        }      
      }
      
    }else{
      posPointFinal=-1
    }
    
    
    length(pp$y)
    
    
    if(posPointFinal==-1){
      return(-1)
    } else{
      return(pp)
    } 
    
    
  }
  
} 



#' @export
getMaximumCurvaturePoints = function(bmdValues,accumulated,limit_curve_up, startPos=1,choiceSmoothing="new",dfreedom=20,scam.m=2){
  
  elbows = NA
  
  originalBmdValues = bmdValues
  originalAccumulated = accumulated
  #plot(originalBmdValues,originalAccumulated,log="x")
  
  # We need to know the length after limiting the curve by the provided threshold
  if(length(bmdValues[bmdValues<=limit_curve_up])<5){
    return(list(NA,NA))
  }
  else{
    pp=getSmoothedCurve(bmdValues,accumulated,limit_curve_up,startPos,choiceSmoothing,dfreedom,scam.m)
    
    #plot(bmdValues,accumulated,log="x")
    #points(10^pp$x,pp$y,col="yellow")
    #plot(pp)
    
    #if pp is equal to -1, we were not able to smooth the line
    if(length(pp)==1){
      elbows_output=list(NA,NA)
    }else{
      posPointFinal = length(pp$x)
      
      bmdValues=10^pp$x[1:posPointFinal]
      accumulated=pp$y[1:posPointFinal]
      
      if(length(bmdValues)<=5){
        elbows_output=list(bmdValues[length(bmdValues)],accumulated[length(bmdValues)])	
      }else{
        elbows_output=findElbow(bmdValues,accumulated,Log=TRUE,S=1)	
        #elbows = as.numeric(elbows)
        
      }
    }
    
    
    return(elbows_output)
    
  }
  
}




#' @export
runPODAccMethod = function(list_bmd_values,choiceSmoothing="flexible",antimode=FALSE){
  
  if(choiceSmoothing=="flexible"){
    
    if(length(list_bmd_values)>1000){
      choiceSmoothing="pp"
    }else{
      choiceSmoothing="new"
    }
    
    
  }
  
  
  # Get accumulation plot values
  output = generateAccumulationValuesFromListBMDValues(list_bmd_values)
  bmd_Values_acc = output[,1]
  accumulated = output[,2]
  #plot(bmd_Values_acc,accumulated,log="x")
  
  # Calculate method for the whole curve first - Needs package scam
  output_elbows = getMaximumCurvaturePoints(bmd_Values_acc,accumulated,bmd_Values_acc[length(bmd_Values_acc)],choiceSmoothing=choiceSmoothing)
  elbows_all_curve = output_elbows[[1]]
  accumulated_all_curve = output_elbows[[2]]
  
  if(antimode==TRUE){
    aux = returnAntimode(list_bmd_values,printHist=FALSE)
    firstmode= aux[1]
    antimode = aux[2]
    
    if(!is.na(antimode)){
      output_elbows = getMaximumCurvaturePoints(bmd_Values_acc,accumulated,10^antimode,choiceSmoothing=choiceSmoothing)
      elbows_anti_mode = output_elbows[[1]]
      accumulated_anti_mode = output_elbows[[2]]
      
      all_elbows = c(elbows_all_curve,elbows_anti_mode)
      accumulated_elbows = c(accumulated_all_curve,accumulated_anti_mode)
      
      mask = is.na(all_elbows)
      
      all_elbows = all_elbows[!mask]
      accumulated_elbows = accumulated_elbows[!mask]
      
    }else{
      
      all_elbows = elbows_all_curve
      accumulated_elbows =accumulated_all_curve
      
    }
    
    if(!is.na(all_elbows[1])){
      # Print plot - choose closest to first mode
      # difffirstmode = abs(all_elbows-10^firstmode)
      elbow = all_elbows[which.min(all_elbows)]
    }else{
      elbow = NA
    }
    
    
  }else{
    
    all_elbows = elbows_all_curve
    accumulated_elbows = accumulated_all_curve
    firstmode=NA
    antimode=NA
    
    
    elbow = all_elbows[which.min(all_elbows)]
    
    # Calculate the absolute differences from 25
    #differences <- abs(accumulated_elbows - 25)
    
    # Find the index of the minimum difference
    #closest_index <- which.min(differences)
    
    #elbow = all_elbows[closest_index]
    if(length(elbow)==0){
      elbow=NA
    }
    
  }
  
  
  
  
  if((is.na(elbow)&(length(bmd_Values_acc)>=5))){
    elbow = median(bmdValues[1:2])
    all_elbows = elbow
  }
  
  
  results = list(podacc=elbow,candidates_podacc=all_elbows,first_mode= 10^firstmode,first_antimode=10^antimode)
  
  
  return(results)    
}



#' @export
runDistMethods = function(distBMDValues,choiceSmoothing="flexible",antimode=FALSE){
  
  
  results_gene_level = runPODAccMethod(distBMDValues)
  
  fifth_perc = as.numeric(quantile(distBMDValues,0.05))
  tenth_perc = as.numeric(quantile(distBMDValues,0.10))
  twenty_five_pos = distBMDValues[25]
  
  
  
  firstmode = results_gene_level$first_mode
  antimode = results_gene_level$first_antimode
  elbow = results_gene_level$podacc
  all_elbows = results_gene_level$candidates_podacc
  
  
  results_gene_level$fifth_perc = fifth_perc
  results_gene_level$twenty_five_pos = twenty_five_pos
  results_gene_level$tenth_perc = tenth_perc
  
  return(results_gene_level)
  
}


#' @export
plotPODAccResults = function(list_bmd_values,results,titleplot="Accumulation Plot Results",xlab_text = "Dose (mg/kg/day)",ylab_text = "Accumulation",legend_pos="right",legend_rel_pos="outside"){
  
  # Get accumulation plot values
  output = generateAccumulationValuesFromListBMDValues(list_bmd_values)
  bmd_Values_acc = output[,1]
  accumulated = output[,2]
  
  
  if(legend_rel_pos=="outside"){
    par(mar = c(5,4,1.4,8))
  }
  
  plot(bmd_Values_acc,accumulated,log="x",main=titleplot,xlab=xlab_text,ylab=ylab_text)  
  
  
  
  
  
  
  # Get accumulation plot values
  output = generateAccumulationValuesFromListBMDValues(list_bmd_values)
  bmd_Values_acc = output[,1]
  accumulated = output[,2]
  
  abline(v=(results$first_mode),col="red", lwd=3)
  abline(v=(results$first_antimode),col="blue", lwd=3)
  
  for(e in 1:length(results$candidates_podacc)){
    abline(v=results$candidates_podacc[e],col="green", lwd=3,lty=2)
  }
  abline(v=results$podacc,col="green", lwd=3)
  
  if(legend_rel_pos=="outside"){
    opar = par(fig=c(0,1,0,1),oma=c(0,0,0,0),mar=c(0,0,0,0),new=TRUE)
    on.exit(par(opar))
    plot(0,0,type='n',bty='n',xaxt='n',yaxt='n')
    
    
    
  }
  
  legend(x=legend_pos,legend=c("First Mode","Antimode",expression("POD"[Accum])),fill=c("red","blue","green"))
  
  
}

#' @export
plotDistMethodResults = function(list_bmd_values,results,titleplot="Accumulation Plot Results",xlab_text = "Dose (mg/kg/day)",ylab_text = "Accumulation",legend_pos="right",legend_rel_pos="outside",GOBP=NA){
  
  # Get accumulation plot values
  output = generateAccumulationValuesFromListBMDValues(list_bmd_values)
  bmd_Values_acc = output[,1]
  accumulated = output[,2]
  
  
  if(legend_rel_pos=="outside"){
    par(mar = c(5,4,1.4,8))
  }
  
  plot(bmd_Values_acc,accumulated,log="x",main=titleplot,xlab=xlab_text,ylab=ylab_text)  
  
  abline(v=(results$first_mode),col="red", lwd=3)
  abline(v=(results$first_antimode),col="blue", lwd=3)
  
  abline(v=results$fifth_perc,col="pink", lwd=3)
  abline(v=results$tenth_perc,col="brown", lwd=3)
  abline(v=results$twenty_five_pos,col="yellow", lwd=3)
  abline(v=GOBP,col="purple", lwd=3)
  
  for(e in 1:length(results$candidates_podacc)){
    abline(v=results$candidates_podacc[e],col="green", lwd=3,lty=2)
  }
  abline(v=results$podacc,col="green", lwd=3)
  
  if(legend_rel_pos=="outside"){
    opar = par(fig=c(0,1,0,1),oma=c(0,0,0,0),mar=c(0,0,0,0),new=TRUE)
    on.exit(par(opar))
    plot(0,0,type='n',bty='n',xaxt='n',yaxt='n')
  }
  

  text_legend = c(expression("POD"[Accum]),"First Mode","Antimode","5Perc","10Perc","Pos25","GOBP")
  colors_fill = c("green","red","blue","pink","brown","yellow","purple")
  values = c(results$podacc,results$first_mode,results$first_antimode,results$fifth_perc,results$tenth_perc,results$twenty_five_pos,GOBP)
  
  colors_fill=colors_fill[!is.na(values)]
  text_legend = text_legend[!is.na(values)]
  
  legend(x="right",legend=text_legend,fill=colors_fill)
  
  
  
}


#' @export
hello_test = function(){
  print("Hello world")
}

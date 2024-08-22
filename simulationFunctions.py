import numpy as np




def zipfWeights(k,s):
	weights = list(range(1,k+1))
	weights = [pow(j,-s) for j in weights]
	total = sum(weights)
	weights = [j/total for j in weights]
	return (np.array(weights),total)

def makeSampleCounts(n,k,truth):
	categories = list(range(1,k+1))
	size = n
	newSamples = np.random.choice(categories,size=size,p=truth)
	counts = np.bincount(newSamples)
	return (counts/n)

def computeEPEerror(sample,truth,epeNormalizer,useEPEnormalizer):
	L = len(sample)
	topProbs = truth[:(L-1)]
	topCounts = sample[1:]
	topError = np.sum(np.abs(topProbs-topCounts))
	
	bottomProbs = truth[L:]
	totalError = topError+np.sum(bottomProbs)
	
	if useEPEnormalizer:
		totalError = totalError/epeNormalizer
	return (totalError)

def makeEPEvector(sample,truth):
	L = len(sample)
	topCounts = sample[1:]
	bottomCounts = np.zeros(len(truth)-(L-1))
	return (np.concatenate([topCounts,bottomCounts]))
	

def computeSSerror(sample,truth):
	L = len(sample)
	countsDecreasingIndex = (-sample).argsort()
	result = 0
	j = 0
	for ind in countsDecreasingIndex:
		if ind != 0:
			est = truth[j]
			trueVal = truth[ind-1]
			result = result+abs(est-trueVal)
			j=j+1
	leftOverZeroCountCategories = np.setdiff1d(list(range(1,len(truth)+1)),countsDecreasingIndex)
	for ind in leftOverZeroCountCategories:
		newError = abs(truth[j]-truth[ind-1])
		result = result+newError
		j=j+1
	return (result)

def makeSSvector(sample,truth):
	result = np.ones(len(truth))
	L = len(sample)
	countsDecreasingIndex = (-sample).argsort()
	j = 0
	for ind in countsDecreasingIndex:
		if ind != 0:
			result[ind-1] = truth[j]
			j=j+1
	leftOverZeroCountCategories = np.setdiff1d(list(range(1,len(truth)+1)),countsDecreasingIndex)
	for ind in leftOverZeroCountCategories:
		if ind != 0:
			result[ind-1] = truth[j]
			j=j+1
	return (result)

	

def makeTSSvector(sample,truth,T):
	result = makeEPEvector(sample,truth)
	ssVector = makeSSvector(sample,truth)
	#First identify some set of T indices such that each index has a count that is at least the T^{th} largest count
	sample = sample[1:]
	countsDecreasingIndex = (-sample).argsort()
	ssTotal = len(countsDecreasingIndex) #ssTotal is both the size of a collection of indices we can use for Sort and Snap and also 1+ the index of the lowest indexed category where there was an observation
	#Now we grab the categories associated with the top T counts
	if ssTotal >= T:
		ssIndices = countsDecreasingIndex[:T]
	else:
		#remaining zero categories
		remIndices = list(range(ssTotal,(len(truth))))
		chosenSSindices = np.array(remIndices[:(T-ssTotal)])
		ssIndices = np.concatenate([countsDecreasingIndex,chosenSSindices])
	result[ssIndices] = ssVector[ssIndices]
	return result

def computeTSSerror(sample,truth,T):
	return (np.sum(np.abs(makeTSSvector(sample,truth,T)-truth)))
		
		
def testMakeVec(n,k,s,T):
	truth = zipfWeights(k,s)
	sample = makeSampleCounts(n,k,truth[0])
	epeVector = makeEPEvector(sample,truth[0])
	ssVector = makeSSvector(sample,truth[0])
	TssVector = makeTSSvector(sample,truth[0],T)
	print("The truth is ")
	print(truth)
	print("The sample is ")
	print(sample)
	print("The epe vector is ")
	print(epeVector)
	print("The ss vector is ")
	print(ssVector)
	print("The Tss vector is ")
	print(TssVector)
		
#testMakeVec(100,15,3,7)
	

#def computeTSSerror(sample,truth,T):
#	#Now use Sort and Snap up to the truncation point
#	L = len(sample)
#	countsDecreasingIndex = (-sample).argsort()
#	result = 0
#	j = 0
#	for ind in countsDecreasingIndex:
#		if ind != 0:
#			est = truth[j]
#			trueVal = truth[ind-1]
#			result = result+abs(est-trueVal)
#			j=j+1
#		if j == T:
#			break
#	#At this point, if j < T, we continue with Sort and Snap. If j = T, we start using empirical proportions at this time
#	leftOverZeroCountCategories = np.setdiff1d(list(range(1,len(truth)+1)),countsDecreasingIndex)
#	if j < T:
#		#in this case, it must be that j = L. First we continue with Sort and Snap until reaching T
#		i = 0
#		for ind in leftOverZeroCountCategories:
#			newError = abs(truth[j]-truth[ind-1])
#			result = result+newError
#			j=j+1
#			i=i+1
#			if j == T:
#				break
#		#Since j > L, all counts are zero.
#		leftOverZeroCountCategories = leftOverZeroCountCategories[i:]
#		result = result+np.sum(truth[leftOverZeroCountCategories-1])
#	if j == T:
#		#in this case, T <= L and we need to apply empirical proportions from here on out
#		if T < L:
#			remainingNonZeroCategories = countsDecreasingIndex[T:]
#			for ind in remainingNonZeroCategories:
#				if ind != 0:
#					est = sample[ind]
#					trueVal = truth[ind-1]
#					result = result+abs(est-trueVal)
#					j=j+1
#	return (result)
	
		

	

def computeError(sample,truth,method,T = 0,epeNormalizer = 0,useEPEnormalizer = False):
	if method == "EPE":
		return computeEPEerror(sample,truth,epeNormalizer,useEPEnormalizer)
	if method == "SS":
		return computeSSerror(sample,truth)
	if method == "TSS":
		return computeTSSerror(sample,truth,T)

def createNs(lowerT,upperT,epsTol,s):
	results = []
	for j in range(lowerT,upperT):
		newN = pow(j,1/(1/(s+2)-epsTol))-1
		results.append(int(newN))
	return results
		

def sim2(dir):
	M = 300
	beta = 1
	ss = [3,5]
	epsTol =.001 #this is for TSS only. Where T is set to n^(1/(s+2)-eps)
	estimators = ["EPE","TSS","SS"]
	eff = open(dir,"w")
	eff.write("n,beta,s,type,mcVal,M\n")
	for s in ss:
		#create Ns
		#Ns = createNs(3,9,epsTol,s) #what is used for s=5
		#if s == 3:
		#	Ns = createNs(7,19,epsTol,3) #what is used for s=3
		Ns = np.exp(list(range(5,15)))
		Ns = Ns.astype(int)
		print("The values of N are ")
		print(Ns)
		for n in Ns:
			for method in estimators:
				print("Now working on n = "+str(n)+":estimator"+method+":s"+str(s))
				k = int(np.floor(pow(n,beta)))
				(truth,total) = zipfWeights(k,s)
				epeNormalizer = zipfWeights(k,s/2)[1]
				#if method is TSS then we compute the truncation point now
				if method == "TSS":
					T = int(np.floor(pow(n,1/(s+2)-epsTol)))
				for m in range(M):
					newSample = makeSampleCounts(n,k,truth)
					if method == "TSS":
						newError = computeError(newSample,truth,method,T,0,False)
					elif method == "EPE":
						newError = computeError(newSample,truth,method,0,epeNormalizer,False)
					elif method == "SS":
						newError = computeError(newSample,truth,method)
					writeString = str(n)+","+str(beta)+","+str(s)+","+str(method)+","+str(newError)+","+str(M)+"\n"
					eff.write(writeString)
	eff.close()

def sim3(dir):
	M = 300
	beta = 1/(1.5)
	ss = [1.5]
	epsTol =.001 #this is for TSS only. Where T is set to n^(1/(s+2)-eps)
	estimators = ["EPE","TSS","SS"]
	eff = open(dir,"w")
	eff.write("n,beta,s,type,mcVal,M\n")
	for s in ss:
		#create Ns
		#Ns = createNs(20,30,epsTol,s)
		Ns = np.exp(list(range(5,15)))
		Ns = Ns.astype(int)
		print("The values of N are ")
		print(Ns)
		for n in Ns:
			for method in estimators:
				print("Now working on n = "+str(n)+":estimator"+method+":s"+str(s))
				k = int(np.floor(pow(n,beta)))
				(truth,total) = zipfWeights(k,s)
				epeNormalizer = zipfWeights(k,s/2)[1]
				#if method is TSS then we compute the truncation point now
				if method == "TSS":
					T = int(np.floor(pow(n,1/(s+2)-epsTol)))
				for m in range(M):
					newSample = makeSampleCounts(n,k,truth)
					if method == "TSS":
						newError = computeError(newSample,truth,method,T,0,False)
					elif method == "EPE":
						newError = computeError(newSample,truth,method,0,epeNormalizer,False)
					elif method == "SS":
						newError = computeError(newSample,truth,method)
					writeString = str(n)+","+str(beta)+","+str(s)+","+str(method)+","+str(newError)+","+str(M)+"\n"
					eff.write(writeString)
	eff.close()


def sim1(dir):
	M = 300
	ss = [1.05]
	beta = 1/(4.05)
	estimators = ["EPE","SS"]
	eff = open(dir,"w")
	eff.write("n,beta,s,type,mcVal,M\n")
	for s in ss:
		#create Ns
		Ns = np.exp(list(range(5,15)))
		Ns = Ns.astype(int)
		print("The values of N are ")
		print(Ns)
		for n in Ns:
			for method in estimators:
				print("Now working on n = "+str(n)+":estimator"+method+":s"+str(s))
				k = int(np.floor(pow(n,beta)))
				(truth,total) = zipfWeights(k,s)
				epeNormalizer = zipfWeights(k,s/2)[1]
				for m in range(M):
					newSample = makeSampleCounts(n,k,truth)
					if method == "EPE":
						newError = computeError(newSample,truth,method,0,epeNormalizer,False)
					elif method == "SS":
						newError = computeError(newSample,truth,method)
					writeString = str(n)+","+str(beta)+","+str(s)+","+str(method)+","+str(newError)+","+str(M)+"\n"
					eff.write(writeString)
	eff.close()

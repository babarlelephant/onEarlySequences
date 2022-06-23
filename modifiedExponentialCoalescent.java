package beast.evolution.tree.coalescent;

public class modifiedCoalescent extends Coalescent {
	// number of cases by onset from December 8 to January 23
	int [] onsetCurve = new int[] {1,0,0,1,4,0,0,4,1,7,3,2,12,6,11,13,12,13,7,26,11,16,24,16,73,63,65,46,78,64,71,
							125,117,187,168,214,281,296,404,458,572,653,667,1254,1222,1552,1937};
	 							
	double serialIntervalInDays = 5.0;
	double serialIntervalInYear = serialIntervalInDays/365.0;

	int day0 = 32+25;  // January 25 (latest sequence in the dataset)
	double ascertainmentRate = 0.15; // divide by 0.15 the number of cases in the epidemic curve
	// standard deviation of Bernouilli random variable with parameter ascertainmentRate
	double sigma = Math.sqrt(ascertainmentRate*(1-ascertainmentRate)*(1-ascertainmentRate)+ascertainmentRate*ascertainmentRate*(1-ascertainmentRate));

	@Override
	public double calculateLogLikelihood(IntervalList intervals, PopulationFunction popSizeFunction, double threshold) {
		ExponentialGrowth eg = (ExponentialGrowth) popSizeFunction;
		double logL = 0;
		for (int j = 0; j < 33; j++) { // we consider the epidemic curve only From Dec 8 to Jan 9
			int day = j+8;
			int daysInPast = day-day0; // compared to day0 = Jan 25 
			double t = -daysInPast/365.0; // in year
			double effectivePopSize = eg.getPopSize(t);
			double modeledPopulationSize = effectivePopSize/serialIntervalInYear;
			double modeledCasesPerDay = modeledPopulationSize/serialIntervalInDays;
			double expectancyOfAscertainedCasePerDay = modeledCasesPerDay*ascertainmentRate;
			double normalizedDiff = Math.abs(onsetCurve[j] - expectancyOfAscertainedCasePerDay) / Math.sqrt(modeledCasesPerDay) / sigma;
			// normalizedDiff follows approximately a standard normal law
			logL -= normalizedDiff * 2 / Math.sqrt(2*Math.PI); // 1st order approximation
		}				
		return logL+super.calculateLogLikelihood(intervals,popSizeFunction,threshold);
	}	
}

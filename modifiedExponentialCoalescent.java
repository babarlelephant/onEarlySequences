package beast.evolution.tree.coalescent;

public class modifiedExponentialCoalescent extends Coalescent {
	
  	// number of cases by onset from December 8 to January 23
	// Yes the SARS-CoV-2 onsetCurve from https://www.thelancet.com/journals/laninf/article/PIIS1473-3099(20)30981-6/fulltext#seccestitle150
  	// is hard-coded into our BEAST2 package
	int [] onsetCurve = new int[] {1,0,0,1,4,0,0,4,1,7,3,2,12,6,11,13,12,13,7,26,11,16,24,16,73,63,65,46,78,64,71,
									125,117,187,168,214,281,296,404,458,572,653,667,1254,1222,1552,1937};
	 							
	
	double serialInterval = 5/365.0; // 5 days
	int day0 = 32+25;  // January 25 (latest sequence in the dataset)
	double ascertainmentRate = 0.15; // divide by 0.15 the number of cases in the epidemic curve
	double sigma = 0.001; 

	@Override
	public double calculateLogLikelihood(IntervalList intervals, PopulationFunction popSizeFunction, double threshold) {
		
		ExponentialGrowth eg = (ExponentialGrowth) popSizeFunction;

		// to help understand what is happening when setting a breakpoint
		double realPopSizeOnJan25 = eg.getPopSize(0)/serialInterval;
		double ascertainedCasesOnJan25=realPopSizeOnJan25*0.15;

		double logL = 0;

		for (int j = 0; j < onsetCurve.length; j++) {
			int day = j+8;
			int daysInPast = day0-day; // compared to day0 = Jan 25 
			double t = daysInPast/365.0; // in year
			double effectivePopSize = eg.getPopSize(t);
			double realPopulationSize = effectivePopSize/serialInterval;

			double normalizedDiff = (onsetCurve[j]/ascertainmentRate)/realPopulationSize-1;
			logL -= normalizedDiff*normalizedDiff*sigma*realPopulationSize;  
			// the cdf is a gaussian, the pdf is  2x exp(-x^2), probably doesn't change much
		}
					
		return logL+super.calculateLogLikelihood(intervals,popSizeFunction,threshold);
	}
}

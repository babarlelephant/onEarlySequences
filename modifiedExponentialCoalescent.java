package beast.evolution.tree.coalescent;



public class modifiedExponentialCoalescent extends Coalescent {
	
	// Yes the SARS-CoV-2 onsetCurve from https://www.thelancet.com/journals/laninf/article/PIIS1473-3099(20)30981-6/fulltext#seccestitle150
  // is hard-code into our BEAST2 package :-)
  
	int [] onsetCurve = new int[] {1,0,0,1,4,0,0,4,1,7,3,2,12,6,11,13,12,13,7,26,11,16,24,16,73,63,65,46,78,64,71,
									125,117,187,168,214,281,296,404,458,572,653,667,1254,1222,1552,1937};
	 							// number of cases by onset from December 8 to January 23
	
	double serialInterval = 5/365.0; // 5 days
	
	/*
	 19 sequences
	 int day0 = 27; // December 27
	 */
	/* 256 sequences */
	int day0 = 32+25;  // January 25 (latest sequence in the dataset)
	double ascertainmentRate = 0.15; // divide by 0.15 the number of cases in the epidemic curve
	double weightParam = 0.1; // the variance of the Gaussians utilised for modeling the cases by onset observations given the modeled exponential curve.

	@Override
	public double calculateLogLikelihood(IntervalList intervals, PopulationFunction popSizeFunction, double threshold) {
			ExponentialGrowth eg = (ExponentialGrowth) popSizeFunction;
			double logL = 0;
			
			double realPopSizeOnJan25 = eg.getPopSize(0)/serialInterval;
			double ascertainedCasesOnJan25=realPopSizeOnJan25*0.15;
			
			for (int j = 0; j < onsetCurve.length; j++) {
				int day = j+8;
				int daysInPast = day0-day; // compared to day0 = Jan 25 
				double t = daysInPast/365.0; // in year
				double effectivePopSize = eg.getPopSize(t);
				double realPopulationSize = effectivePopSize/serialInterval;
				
				// too naive, not quite the best modeling
				double normalizedDiff = (onsetCurve[j]/ascertainmentRate)/realPopulationSize-1;
				logL -= normalizedDiff*normalizedDiff * weightParam; 
			}
				
					
			return logL+super.calculateLogLikelihood(intervals,popSizeFunction,threshold);
	}
	


	
}

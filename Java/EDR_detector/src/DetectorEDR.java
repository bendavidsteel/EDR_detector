import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

import Jama.Matrix;

public class DetectorEDR {

	static double mutation_rate = 0.1;
	static int population = 20;
	static int no_genes = 9;
	static int no_parents = 2;
	
	
	public static void main(String[] args) {
		String path = "C:/Projects/Matlab/SCR_labelling/";
		
		System.out.println("Loading CSV data...");
		
		List<List<double[]>> edaCSVs= new ArrayList<>();
		List<List<double[]>> targetCSVs = new ArrayList<>();
		
		List<double[]> all_data = new ArrayList<>();
		List<long[]> all_times = new ArrayList<>();
		
		List<long[]> targets = new ArrayList<>();
		
		for (int j = 1; j <= 3; j++) {
			String csvFile = "eda_data/recording_config" + String.valueOf(j) + "_eda.csv";
			
			edaCSVs.add(getCSV(path + csvFile));
			
			String targetFile = "peak_locations_config" + String.valueOf(j) + ".txt";
			
			targetCSVs.add(getCSV(path + targetFile));
			
			double[] temp_data = new double[edaCSVs.get(j-1).size()];
			long[] temp_times = new long[edaCSVs.get(j-1).size()];
			
			for (int i = 0; i < edaCSVs.get(j-1).size(); i++) {
				temp_times[i] = (long) edaCSVs.get(j-1).get(i)[0];
				temp_data[i] = edaCSVs.get(j-1).get(i)[1];
			}
			
			all_data.add(temp_data);
			all_times.add(temp_times);
			
			long[] temp_targets = new long[targetCSVs.get(j-1).size()];
			
			for (int i = 0; i < targetCSVs.get(j-1).size(); i++) {
				temp_targets[i] = (long) targetCSVs.get(j-1).get(i)[0];
			}
			
			targets.add(temp_targets);
		}
		
		System.out.println("Generating Initial and Running Population...");
		
		Gene[][] chromosomes = generateInitialPopulation();
		
		int[] scores = runPopulation(all_data, all_times, chromosomes, targets);
		
		double[][] parents = new double[no_genes][no_parents];
		
		int generation = 0;
		
		int greatest_algorithm = -100;
		
		System.out.println("Starting Iterating through Generations...");
		
		while (generation < 10000000) {
			
			
			//Selecting parents of next generation
			int[] indices = bubbleSortIndices(scores);
			
			int generational_score = scores[indices[indices.length-1]];
			
			if (generational_score > greatest_algorithm) {
				
				greatest_algorithm = generational_score;
				
				String best_parameters = "";
				
				for (int i = 0; i < no_genes; i++) {
					best_parameters += String.valueOf(chromosomes[i][indices[indices.length-1]].value) + " ";
				}
				
				List<String> lines = Arrays.asList(best_parameters);
				Path file = Paths.get("C:/Users/benda/eclipse-workspace/EDR_detector/algorithm.txt");
				try {
					Files.write(file, lines, Charset.forName("UTF-8"));
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
			
			System.out.println("Generation: " + String.valueOf(generation));
			System.out.println("Score: " + String.valueOf(generational_score));
			
			for (int i = 0; i < no_parents; i++) {
				for (int j = 0; j < no_genes; j++) {
					parents[j][i] = chromosomes[j][indices[population - i - 1]].value;
				}
			}
			
			
			//Breeding new generation
			Random rand = new Random();
			
			for (int i = 0; i < population; i++) {
				for (int j = 0; j < no_genes; j++) {
					if (rand.nextFloat() < mutation_rate) {
						chromosomes[j][i].mutate();
					} else {
						chromosomes[j][i].value = parents[j][rand.nextInt(no_parents)];
					}
				}
			}
			
			scores = runPopulation(all_data, all_times, chromosomes, targets);
			
			generation += 1;
		}
	}
	
	public static List<double[]> getCSV(String file) {

        String csvFile = file;
        String line = "";
        String cvsSplitBy = ";";
        
        List<double[]> dataCSV = new ArrayList<>();

        try (BufferedReader br = new BufferedReader(new FileReader(csvFile))) {

            while ((line = br.readLine()) != null) {

                // use comma as separator
                String[] row = line.split(cvsSplitBy);
                
                double[] dataRow = new double[row.length];
                
                for (int i = 0; i < row.length; i++) {
                	dataRow[i] = Double.valueOf(row[i].replaceAll("\"", ""));
                }
                
                dataCSV.add(dataRow);

            }

        } catch (IOException | NumberFormatException e) {
            e.printStackTrace();
        }
        
        return dataCSV;
    }
	
	public static Gene[][] generateInitialPopulation() {
		
		Gene[][] chromosomes = new Gene[no_genes][population];
		
		for (int i = 0; i < population; i++) {
			//Assigning randomly generated starting values to each gene in chromosome
			
			Gene len_data = new Gene(10, 40);
			chromosomes[0][i] = len_data;
			
			Gene overlap = new Gene(1, 5);
			chromosomes[1][i] = overlap;
			
			Gene offset = new Gene(0, 40);
			chromosomes[2][i] = offset;
			
			Gene start_WT = new Gene(0, 40);
			chromosomes[3][i] = start_WT;
			
			Gene end_WT = new Gene(0, 40);
			chromosomes[4][i] = end_WT;
			
			Gene thres_low = new Gene(0, 1);
			chromosomes[5][i] = thres_low;
			
			Gene theta = new Gene(1, 0.1);
			chromosomes[6][i] = theta;
			
			Gene alpha = new Gene(0, 1);
			chromosomes[7][i] = alpha;
			
			Gene width = new Gene(1, 10);
			chromosomes[8][i] = width;
		}
		
		return chromosomes;
	}
	
	
	public static int[] runPopulation(List<double[]> all_data, List<long[]> all_times, Gene[][] chromosomes, List<long[]> targets) {
		
		int[] scores = new int[population];
		
		for (int i = 0; i < population; i++) {
			for (int j = 0; j < all_data.size(); j++) {
			
				List<SCR> peaks = realTimeFinder(all_data.get(j), all_times.get(j), chromosomes[0][i].value,
																		  			chromosomes[1][i].value,
																		  			chromosomes[2][i].value,
																		  			chromosomes[3][i].value,
																		  			chromosomes[4][i].value,
																		  			chromosomes[5][i].value,
																		  			chromosomes[6][i].value,
																		  			chromosomes[7][i].value,
																		  			chromosomes[8][i].value);
				
				scores[i] += fitness(targets.get(j), peaks, chromosomes[0][i].value, chromosomes[1][i].value);
			}
			
		}	
		
		return scores;
	}
	
	public static List<SCR> realTimeFinder(double[] all_data, long[] all_times, double len_data_d, 
																					   double overlap_d,
																					   double offset_d, 
																					   double start_WT_d,
																					   double end_WT_d,
																					   double thres_low,
																					   double theta,
																					   double alpha,
																					   double width_d) {
		
		int len_data = (int) Math.round(len_data_d);
		int overlap = (int) Math.round(overlap_d);
		int offset = (int) Math.round(offset_d);
		int start_WT = (int) Math.round(start_WT_d);
		int end_WT = (int) Math.round(end_WT_d);
		int width = (int) Math.round(width_d);
		
		List<SCR> peaks = new ArrayList<>();
		
		for (int i = (int) len_data; i < all_data.length; i++) {
			if (all_times[i] % (len_data/overlap) == 0) {
				long[] times = Arrays.copyOfRange(all_times, i - len_data, i);
				double[] data = Arrays.copyOfRange(all_data, i - len_data, i);
				
				double[] median_data = data;
				if (width > 1) {
					median_data = medianFilter(data, width);
				}
				
				//double[] poly = polynomialRegression(data, range, poly_degree);
				
				//double[] poly_deriv = calculateDerivative(poly);
				
				//double[] data_deriv = evaluatePolynomial(times, poly_deriv);
			
				
				List<SCR> SCRs = findPeaks2(toList(median_data), toList(times), offset,
																					start_WT, 
																					end_WT,
																					thres_low,
																					theta,
																					alpha);
				
				for (SCR scr : SCRs) {
					boolean newSCR = true;
					
					for (SCR peak : peaks) {
						if (peak.peak_time == scr.peak_time) {
							newSCR = false;
						}
					}
					
					if (newSCR) {
						peaks.add(scr);
					}
				}
			}
		}
		return peaks;
	}
	
	public static int fitness(long[] targets, List<SCR> peaks, double len_data, double overlap) {
		int score = 0;
		
		int[] matches = new int[peaks.size()];
	    int[] target_matches = new int[targets.length];
	    
	    if (peaks.size() == 0) {
	        score = -100;
	    } else {
	        for (int i = 0; i < peaks.size(); i++) {
	            matches[i] = 0;
	            for (int j = 0; j < targets.length; j++) {
	            	
	                int diff = (int) Math.abs(targets[j] - peaks.get(i).peak_time);
	                
	                if (diff < 100) {
	                    score += 10;
	                    matches[i] = j;
	                    target_matches[j] = 1;
	                } else if (diff < 200) {
	                    score += 5;
	                    matches[i] = j;
	                    target_matches[j] = 1;
	                }
	            }
	        

	            if (matches[i] == 0) {
	                score -= 20;
	            }
	        }

	        int[] hist = new int[targets.length];
	        
	        for (int i = 0; i < matches.length; i++) {
	            if (matches[i] != 0) {
	                hist[matches[i]] = hist[matches[i]] + 1;
	            }
	        }
	        for (int i = 0; i < hist.length; i++) {
	            if (hist[i] > 1) {
	                score -= 10;
	            }
	        }
	        
	        score = score + 5*sum(target_matches);
	        
	        score -= (int) 0.1*len_data;
	        score -= (int) 1*overlap/len_data;
		}
		
		return score;
	}
	
	
	public static List<SCR> findPeaks(List<Double> data, double[] data_deriv, List<Long> times, int offset, 
																			   				   		   int start_WT, 
																			   				   		   int end_WT, 
																			   				   		   double thres_low, 
																			   				   		   double theta) {

        /*
        This function finds the peaks of an EDA signal and returns basic properties.
		Also, peak_end is assumed to be no later than the start of the next peak. (Is this okay??)

		********* INPUTS **********
		data:        DataFrame with EDA as one of the columns and indexed by a datetimeIndex
		offset:      the number of rising samples and falling samples after a peak needed to be counted as a peak
		start_WT:    maximum number of seconds before the apex of a peak that is the "start" of the peak
		end_WT:      maximum number of seconds after the apex of a peak that is the "rec.t/2" of the peak, 50% of amp
		thres:       the minimum uS change required to register as a peak, defaults as 0 (i.e. all peaks count)
		sampleRate:  number of samples per second, default=8
         */


		int sampleRate = 15;
	

        int[] peaks = new int[data.size()];
        Arrays.fill(peaks, 0);

        long[] peak_times = new long[data.size()];
        Arrays.fill(peak_times, 0);

        double[] peak_eda = new double[data.size()];
        Arrays.fill(peak_eda, 0);

        int[] peak_sign = new int[data_deriv.length];

        for (int i = 0; i < data_deriv.length; i++) {
            if (data_deriv[i] > 0) {
                peak_sign[i] = 1;
            } else if (data_deriv[i] < 0) {
                peak_sign[i] = -1;
            } else {
                peak_sign[i] = 0;
            }
        }

        for (int i = (int) offset; i < (data_deriv.length - offset - 1); i++) {
            if ((peak_sign[i] == 1) && (peak_sign[i + 1] < 1)) {

                peaks[i] = 1;
                peak_times[i] = times.get(i);
                peak_eda[i] = data.get(i);

                for (int j = 1; j < offset; j++) {
                    if ((i - j > 0) && (i + j > peak_sign.length)) {
                        if ((peak_sign[i - j] < 1) || (peak_sign[i + j] > -1)) {
                            peaks[i] = 0;
                            break;
                        }
                    }
                }
            }
        }


        int[] peak_start = new int[data.size()];
        Arrays.fill(peak_start, 0);

        long[] peak_start_times = new long[data.size()];
        Arrays.fill(peak_start_times, 0);

        double[] max_deriv = new double[data.size()];
        Arrays.fill(max_deriv, 0);

        long[] rise_time = new long[data.size()];
        Arrays.fill(rise_time, 0);


        //Finding start of peaks
        for (int i = 0; i < peaks.length - 1; i++) {
            if (peaks[i] == 1) {
                int temp_start = 0;
                if ((i - sampleRate) > 0) {
                    temp_start = i - sampleRate;
                }

                max_deriv[i] = max(Arrays.copyOfRange(data_deriv, temp_start, i));
                double start_deriv = theta * max_deriv[i];

                boolean found = false;
                int find_start = i;

                //Has to peak within start_WT seconds
                while ((!found) && (find_start > (i - (start_WT * sampleRate))) && (find_start >= 0)) {
                    if (data_deriv[find_start] < start_deriv) {
                        found = true;
                        peak_start[find_start] = 1;
                        peak_start_times[i] = times.get(find_start);
                        rise_time[i] = times.get(i) - peak_start_times[i];
                    }

                    find_start -= 1;
                }

                //If we didn't find a start
                if (!found) {
                    if ((i - (start_WT * sampleRate)) >= 0) {
                        peak_start[i - ((int) start_WT * sampleRate)] = 1;
                        peak_start_times[i] = times.get(i - ((int) start_WT * sampleRate));
                        rise_time[i] = start_WT;
                    } else {
                        peak_start[0] = 1;
                        peak_start_times[i] = times.get(0);
                        rise_time[i] = start_WT;
                    }
                }

                if ((thres_low > 0) && (data.get(i) - data.get(times.indexOf(peak_start_times[i])) < thres_low)) {
                    peaks[i] = 0;
                    peak_start[i] = 0;
                    peak_start_times[i] = 0;
                    max_deriv[i] = 0;
                    rise_time[i] = 0;
                    peak_times[i] = 0;
                    peak_eda[i] = 0;
                }
            }
        }

        int[] peak_end = new int[data.size()];
        Arrays.fill(peak_end, 0);

        long[] peak_end_times = new long[data.size()];
        Arrays.fill(peak_end_times, 0);

        double[] amplitude = new double[data.size()];
        Arrays.fill(amplitude, 0);

        float[] decay_time = new float[data.size()];
        Arrays.fill(decay_time, 0);

        long[] half_rise = new long[data.size()];
        Arrays.fill(half_rise, 0);

        float[] SCR_width = new float[data.size()];
        Arrays.fill(SCR_width, 0);

        //Finding the end of the peak, amplitude of peak
        for (int i = 0; i < peaks.length - 1; i++) {
            if (peaks[i] == 1) {
                double peak_amp = data.get(i);
                double start_amp = data.get(times.indexOf(peak_start_times[i]));
                amplitude[i] = peak_amp - start_amp;

                double half_amp = (amplitude[i] * 0.5) + start_amp;

                boolean found = false;
                int find_end = i;

                //Has to decay within end_WT seconds
                while ((!found) && (find_end < (i + (end_WT * sampleRate))) && (find_end < (peaks.length))) {
                    if (data.get(find_end) < half_amp) {
                        found = true;
                        peak_end[find_end] = 1;
                        peak_end_times[i] = Math.round(times.get(find_end));
                        decay_time[i] = (float) (peak_end_times[i] - times.get(i));

                        //Find width
                        int find_rise = i;
                        boolean found_rise = false;

                        while (!found_rise) {
                            if (data.get(find_rise) < half_amp) {
                                found_rise = true;
                                half_rise[i] = Math.round(times.get(find_rise));
                                SCR_width[i] = (float) (peak_end_times[i] - times.get(find_rise));
                            }

                            find_rise -= 1;
                        }
                    } else if (peak_start[find_end] == 1) {
                        found = true;
                        peak_end[find_end] = 1;
                        peak_end_times[i] = Math.round(times.get(find_end));
                    }

                    find_end -= 1;
                }

                if (!found) {
                    int min_index = minIndex(data.subList(i, i + ((int) end_WT * sampleRate)));
                    peak_end[i + min_index] = 1;
                    peak_end_times[i] = Math.round(times.get(i + min_index));
                }
            }
        }

        for (int i = 0; i < max_deriv.length; i++) {
            max_deriv[i] *= sampleRate;
        }

        List<SCR> result = new ArrayList<>();

        for (int i = 0; i < data.size(); i++) {
            if (peaks[i] == 1) {
                SCR scr = new SCR(peak_times[i], peak_eda[i], rise_time[i], max_deriv[i], amplitude[i]);

                result.add(scr);
            }
        }
        return result;
    }
	
	public static List<SCR> findPeaks2(List<Double> data, List<Long> times, int offset, 
	   		   																	 int start_WT, 
	   		   																	 int end_WT, 
	   		   																	 double thres_low, 
	   		   																	 double theta,
	   		   																	 double alpha) {

		/*
		This function finds the peaks of an EDA signal and returns basic properties.
		Also, peak_end is assumed to be no later than the start of the next peak. (Is this okay??)
		
		********* INPUTS **********
		data:        DataFrame with EDA as one of the columns and indexed by a datetimeIndex
		offset:      the number of rising samples and falling samples after a peak needed to be counted as a peak
		start_WT:    maximum number of seconds before the apex of a peak that is the "start" of the peak
		end_WT:      maximum number of seconds after the apex of a peak that is the "rec.t/2" of the peak, 50% of amp
		thres:       the minimum uS change required to register as a peak, defaults as 0 (i.e. all peaks count)
		sampleRate:  number of samples per second, default=8
		*/
		
		
		int sampleRate = 15;
		
		double[] filtered_data = toArrayDouble(data);
		
		for (int i = 1; i < data.size(); i++) {
			filtered_data[i] = (alpha * filtered_data[i-1]) + ((1 - alpha) * filtered_data[i]);
		}
		
		double[] data_deriv = new double[filtered_data.length - 1];
		
		for (int i = 0; i < data_deriv.length; i++) {
			data_deriv[i] = filtered_data[i+1] - filtered_data[i];
		}
		
		int[] peaks = new int[data.size()];
		Arrays.fill(peaks, 0);
		
		long[] peak_times = new long[data.size()];
		Arrays.fill(peak_times, 0);
		
		double[] peak_eda = new double[data.size()];
		Arrays.fill(peak_eda, 0);
		
		int[] peak_sign = new int[data_deriv.length];
		
		for (int i = 0; i < data_deriv.length; i++) {
			if (data_deriv[i] > 0) {
				peak_sign[i] = 1;
			} else if (data_deriv[i] < 0) {
				peak_sign[i] = -1;
			} else {
				peak_sign[i] = 0;
			}
		}
		
		for (int i = (int) offset; i < (data_deriv.length - offset - 1); i++) {
			if ((peak_sign[i] == 1) && (peak_sign[i + 1] < 1)) {
			
				peaks[i] = 1;
				peak_times[i] = times.get(i);
				peak_eda[i] = data.get(i);
				
				for (int j = 1; j < offset; j++) {
					if ((i - j > 0) && (i + j > peak_sign.length)) {
						if ((peak_sign[i - j] < 1) || (peak_sign[i + j] > -1)) {
							peaks[i] = 0;
							break;
						}
					}
				}
			}
		}
		
		
		int[] peak_start = new int[data.size()];
		Arrays.fill(peak_start, 0);
		
		long[] peak_start_times = new long[data.size()];
		Arrays.fill(peak_start_times, 0);
		
		double[] max_deriv = new double[data.size()];
		Arrays.fill(max_deriv, 0);
		
		long[] rise_time = new long[data.size()];
		Arrays.fill(rise_time, 0);
		
		
		//Finding start of peaks
		for (int i = 0; i < peaks.length - 1; i++) {
			if (peaks[i] == 1) {
				int temp_start = 0;
				if ((i - sampleRate) > 0) {
					temp_start = i - sampleRate;
				}
				
				max_deriv[i] = max(Arrays.copyOfRange(data_deriv, temp_start, i));
				double start_deriv = theta * max_deriv[i];
				
				boolean found = false;
				int find_start = i;
				
				//Has to peak within start_WT seconds
				while ((!found) && (find_start > (i - (start_WT * sampleRate))) && (find_start >= 0)) {
					if (data_deriv[find_start] < start_deriv) {
						found = true;
						peak_start[find_start] = 1;
						peak_start_times[i] = times.get(find_start);
						rise_time[i] = times.get(i) - peak_start_times[i];
					}
				
					find_start -= 1;
				}
				
				//If we didn't find a start
				if (!found) {
					if ((i - (start_WT * sampleRate)) >= 0) {
						peak_start[i - ((int) start_WT * sampleRate)] = 1;
						peak_start_times[i] = times.get(i - ((int) start_WT * sampleRate));
						rise_time[i] = start_WT;
					} else {
						peak_start[0] = 1;
						peak_start_times[i] = times.get(0);
						rise_time[i] = start_WT;
					}
				}
				
				if ((thres_low > 0) && (data.get(i) - data.get(times.indexOf(peak_start_times[i])) < thres_low)) {
					peaks[i] = 0;
					peak_start[i] = 0;
					peak_start_times[i] = 0;
					max_deriv[i] = 0;
					rise_time[i] = 0;
					peak_times[i] = 0;
					peak_eda[i] = 0;
				}
			}
		}
		
		int[] peak_end = new int[data.size()];
		Arrays.fill(peak_end, 0);
		
		long[] peak_end_times = new long[data.size()];
		Arrays.fill(peak_end_times, 0);
		
		double[] amplitude = new double[data.size()];
		Arrays.fill(amplitude, 0);
		
		float[] decay_time = new float[data.size()];
		Arrays.fill(decay_time, 0);
		
		long[] half_rise = new long[data.size()];
		Arrays.fill(half_rise, 0);
		
		float[] SCR_width = new float[data.size()];
		Arrays.fill(SCR_width, 0);
		
		//Finding the end of the peak, amplitude of peak
		for (int i = 0; i < peaks.length - 1; i++) {
			if (peaks[i] == 1) {
				double peak_amp = data.get(i);
				double start_amp = data.get(times.indexOf(peak_start_times[i]));
				amplitude[i] = peak_amp - start_amp;
				
				double half_amp = (amplitude[i] * 0.5) + start_amp;
				
				boolean found = false;
				int find_end = i;
				
				//Has to decay within end_WT seconds
				while ((!found) && (find_end < (i + (end_WT * sampleRate))) && (find_end < (peaks.length))) {
					if (data.get(find_end) < half_amp) {
						found = true;
						peak_end[find_end] = 1;
						peak_end_times[i] = Math.round(times.get(find_end));
						decay_time[i] = (float) (peak_end_times[i] - times.get(i));
						
						//Find width
						int find_rise = i;
						boolean found_rise = false;
						
						while (!found_rise) {
							if (data.get(find_rise) < half_amp) {
								found_rise = true;
								half_rise[i] = Math.round(times.get(find_rise));
								SCR_width[i] = (float) (peak_end_times[i] - times.get(find_rise));
							}
						
						find_rise -= 1;
						}
					} else if (peak_start[find_end] == 1) {
						found = true;
						peak_end[find_end] = 1;
						peak_end_times[i] = Math.round(times.get(find_end));
					}
					
					find_end -= 1;
				}
				
				if (!found) {
					int min_index = minIndex(data.subList(i, i + ((int) end_WT * sampleRate)));
					peak_end[i + min_index] = 1;
					peak_end_times[i] = Math.round(times.get(i + min_index));
				}
			}
		}
		
		for (int i = 0; i < max_deriv.length; i++) {
			max_deriv[i] *= sampleRate;
		}
		
		List<SCR> result = new ArrayList<>();
		
		for (int i = 0; i < data.size(); i++) {
			if (peaks[i] == 1) {
				SCR scr = new SCR(peak_times[i], peak_eda[i], rise_time[i], max_deriv[i], amplitude[i]);
	
	            result.add(scr);
			}
		}
		return result;
	}
	
	
	public static int[] findResponses(ArrayList<Float> lastX2Deriv, ArrayList<Float> lastXData, ArrayList<Long> lastXTime, int polyDegree) {

        int flag[] = null;

        int maxIndex = 0;
        int minIndex = 0;

        for (int i = 0; i < lastX2Deriv.size(); i++) {
            if (lastX2Deriv.get(i) > lastX2Deriv.get(maxIndex)) {
                maxIndex = i;
            }
            if (lastX2Deriv.get(i) < lastX2Deriv.get(minIndex)) {
                minIndex = i;
            }
        }

        if ((minIndex > maxIndex) && (lastX2Deriv.get(maxIndex) > 0) && (lastX2Deriv.get(minIndex) < 0)) {
            int zeroCrossing = maxIndex;

            for (int i = 0; i < minIndex; i++) {
                if (Math.abs(lastX2Deriv.get(i)) < Math.abs(lastX2Deriv.get(zeroCrossing))) {
                    zeroCrossing = i;
                }
            }

            flag = new int[3];
            flag[0] = maxIndex;
            flag[1] = zeroCrossing;
            flag[2] = minIndex;
        }

        return flag;
    }
	
	
	public static double[] polynomialRegression(double[] originalData,double[] timeInMilli, int degree){
        double data[] = originalData;
        int n = data.length;
        double tempTime[] = timeInMilli;
        double x_matrix[][] = new double[degree+1][n];

        //Here we calculate the x matrix
        for(int i = 0 ; i < n ; i++) {
            for (int j = 0; j < degree+1; j++) {
                x_matrix[j][i] = Math.pow(tempTime[i] - tempTime[n-1], j);
            }
        }

        Matrix A = new Matrix(x_matrix);
        Matrix y = new Matrix(data,1);

        //Here we use the matric formula to calculate the coefficients
        Matrix Ainv = A.transpose().inverse();
        Matrix x = Ainv.times(y.transpose());


        Matrix new_y = x.transpose().times(A);
        double coefficient[][] = x.getArray();
        return flattenArray(coefficient);
    }
	
	
	public static double[] calculateDerivative(double coefficients[]) {
        double derivatedCoefficient[] = new double[coefficients.length - 1];
        double invertedCoefficient[] = new double[coefficients.length];
        for (int i = 0; i < coefficients.length; i++) {
            invertedCoefficient[i] = coefficients[coefficients.length - 1 - i];
        }

        for (int i = 0; i < invertedCoefficient.length - 1; i++) {
            derivatedCoefficient[derivatedCoefficient.length - 1 - i] = invertedCoefficient[i] * (invertedCoefficient.length - i - 1);
        }
        return derivatedCoefficient;
    }
	
	public static double[] evaluatePolynomial(long[] x, double[] poly) {
		
		double[] y = new double[x.length];
		
		for (int i = 0; i < x.length; i++) {
			for (int j = 0; j < poly.length; j++) {
				y[i] += poly[j]*Math.pow(x[i],(poly.length - 1 - i));
			}
		}
		
		return y;
	}
	
	public static double[] medianFilter(double[] data, int width) {
		
		double[] median_data = new double[data.length];
		
		for (int i = 0; i < data.length; i++) {
			int start;
			
			if (i - (int)((double)width / 2) >= 0) {
				start = i - (int)((double)width / 2);
			} else {
				start = 0;
			}
			
			int finish;
			
			if (i + (int)((double)width / 2) + 1 < data.length) {
				finish = i + (int)((double)width / 2) + 1;
			} else {
				finish = data.length;
			}
			
			int[] sortedIndices = bubbleSortIndices(Arrays.copyOfRange(data, start, finish));
			
			if (sortedIndices.length % 2 != 0) {
				median_data[i] = data[sortedIndices[(int)((double)sortedIndices.length) / 2]];
			} else {
				median_data[i] = (data[sortedIndices[(int)((double)sortedIndices.length / 2)]] + 
								  data[sortedIndices[(int)((double)sortedIndices.length / 2) - 1]]) / 2;
			}
		}
		
		return median_data;
	}

    /////////////////////////////Useful Array Functions/////////////////////////////////

    private static int minIndex(List<Double> array) {
        int minIndex = 0;
        for (int i = 1; i < array.size(); i++) {
            if (array.get(i) < array.get(minIndex)) {
                minIndex = i;
            }
        }
        return minIndex;
    }
    
    private static double max(double[] array) {
    	double max = -Double.MAX_VALUE;
    	
    	for (int i = 0; i < array.length; i++) {
    		if (array[i] > max) {
    			max = array[i];
    		}
    	}
    	
    	return max;
    }
    
    private static int max(int[] array) {
    	int max = -Integer.MAX_VALUE;
    	
    	for (int i = 0; i < array.length; i++) {
    		if (array[i] > max) {
    			max = array[i];
    		}
    	}
    	
    	return max;
    }
    
    public static double[] ArrFloat2float(List<Double> arrData){
        double[] arrResults = new double[arrData.size()];
        for(int i = 0; i < arrData.size(); ++i) {
            arrResults[i] = arrData.get(i);
        }
        return arrResults;
    }

    /*
    * Helper function to reduce an multidimensional array to
    * a smaller dimension that doesn't contain empty row or column
     */
    public static double[] flattenArray(double data[][]) {
        double newArray[] = new double[data.length*data[0].length];
        for(int i = 0; i < data.length; i++) {
            double[] row = data[i];
            for(int j = 0; j < row.length; j++) {
                double number = data[i][j];
                newArray[i*row.length+j] = number;
            }
        }

        return newArray;
    }
    
    public static List<Double> toList(double[] array) {
    	List<Double> list = new ArrayList<>();
    	
    	for (int i = 0; i < array.length; i++) {
    		list.add(array[i]);
    	}
    	
    	return list;
    }
    
    public static List<Long> toList(long[] array) {
    	List<Long> list = new ArrayList<>();
    	
    	for (int i = 0; i < array.length; i++) {
    		list.add(array[i]);
    	}
    	
    	return list;
    }
    
    public static long[] toArrayLong(List<Double> list) {
    	long[] array = new long[list.size()];
    	
    	for (int i = 0; i < array.length; i++) {
    		array[i] = Math.round(list.get(i));
    	}
    	
    	return array;
    }
    
    public static double[] toArrayDouble(List<Double> list) {
    	double[] array = new double[list.size()];
    	
    	for (int i = 0; i < array.length; i++) {
    		array[i] = list.get(i);
    	}
    	
    	return array;
    }
    
    public static int sum(int[] array) {
    	int sum = 0;
    	
    	for (int i = 0; i < array.length; i++) {
    		sum += array[i];
    	}
    	
    	return sum;
    }
    
    public static int[] bubbleSortIndices(int[] array) {
    	
    	int[] indices = new int[array.length];
    	
    	for (int i = 0; i < array.length; i++) {
    		indices[i] = i;
    	}
    	
        boolean swapped = true;
        int j = 0;
        int tmp;
        while (swapped) {
            swapped = false;
            j++;
            for (int i = 0; i < array.length - j; i++) {
                if (array[i] > array[i + 1]) {
                	
                    tmp = array[i];
                    array[i] = array[i + 1];
                    array[i + 1] = tmp;
                    
                    tmp = indices[i];
                    indices[i] = indices[i + 1];
                    indices[i + 1] = tmp;
                    
                    swapped = true;
                }
            }
        }
        
        return indices;
    }
    
    public static int[] bubbleSortIndices(double[] array) {
    	
    	int[] indices = new int[array.length];
    	
    	for (int i = 0; i < array.length; i++) {
    		indices[i] = i;
    	}
    	
        boolean swapped = true;
        int j = 0;
        
        int tmpI;
        double tmpD;
        
        while (swapped) {
            swapped = false;
            j++;
            for (int i = 0; i < array.length - j; i++) {
                if (array[i] > array[i + 1]) {
                	
                    tmpD = array[i];
                    array[i] = array[i + 1];
                    array[i + 1] = tmpD;
                    
                    tmpI = indices[i];
                    indices[i] = indices[i + 1];
                    indices[i + 1] = tmpI;
                    
                    swapped = true;
                }
            }
        }
        
        return indices;
    }
}

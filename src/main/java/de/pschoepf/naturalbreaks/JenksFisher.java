/*
* Port of Jenks/Fisher breaks originally created in C by Maarten Hilferink.
*
* Copyright (C) {2015}  {Philipp Schoepf}
*
*    This program is free software: you can redistribute it and/or modify
*    it under the terms of the GNU General Public License as published by
*    the Free Software Foundation, either version 3 of the License, or
*    (at your option) any later version.
*
*    This program is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*    GNU General Public License for more details.
*
*    You should have received a copy of the GNU General Public License
*    along with this program.  If not, see <http://www.gnu.org/licenses/>.
**/
package de.pschoepf.naturalbreaks;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;

/**
 *
 * Basic port of original C code from Maarten Hilferink.
 * All credits for this fantastic work go to him.
 *
 *
 *
 * @author Philipp Schöpf
 */
public class JenksFisher {

    private List<ValueCountPair> cumulValues;
    private int numValues = 0;
    private int numBreaks = 0;
    private int bufferSize = 0;
    private double[] previousSSM;
    private double[] currentSSM;
    private int[] classBreaks;
    private int classBreaksIndex = 0;
    private int completedRows = 0;

    /**
     * Constructor that initializes main variables used in fisher calculation of natural breaks.
     *
     * @param vcpc Ordered list of pairs of values to occurrence counts.
     * @param k Number of breaks to find.
     */
    public JenksFisher(List<ValueCountPair> vcpc, int k) {
        this.cumulValues = new ArrayList<>();
        this.numValues = vcpc.size();
        this.numBreaks = k;
        this.bufferSize = (vcpc.size() - (k - 1));
        this.previousSSM = new double[this.bufferSize];
        this.currentSSM = new double[this.bufferSize];
        this.classBreaks = new int[this.bufferSize * (this.numBreaks - 1)];
        this.classBreaksIndex = 0;
        this.completedRows = 0;
        double cwv = 0.0;
        int cw = 0, w = 0;

        ValueCountPair currPair;

        for (int i = 0; i != this.numValues; ++i) {
            currPair = vcpc.get(i);
            assert (i == 0 || currPair.getValue() >= vcpc.get(i - 1).getValue()); // PRECONDITION: the value sequence must be strictly increasing

            w = currPair.getCount();
            assert (w > 0); // PRECONDITION: all weights must be positive

            cw += w;
            assert (cw >= w); // No overflow? No loss of precision?

            cwv += w * currPair.getValue();
            this.cumulValues.add(new ValueCountPair(cwv, cw));
            if (i < this.bufferSize) {
                this.previousSSM[i] = cwv * cwv / cw; // prepare sum of squared means for first class. Last (k-1) values are omitted
            }
        }

    }

    /**
     * Gets sum of weighs for elements with index b..e.
     *
     * @param b index of begin element
     * @param e index of end element
     * @return sum of weights.
     */
    private int getSumOfWeights(int b, int e) {
        assert (b != 0);    // First element always belongs to class 0, thus queries should never include it.
        assert (b <= e);
        assert (e < this.numValues);

        int res = this.cumulValues.get(e).getCount();
        res -= this.cumulValues.get(b - 1).getCount();
        return res;
    }

    /**
     * Gets sum of weighed values for elements with index b..e
     *
     * @param b index of begin element
     * @param e index of end element
     * @return the cumul. sum of the values*weight
     */
    private double getSumOfWeightedValues(int b, int e) {
        assert (b != 0);
        assert (b <= e);
        assert (e < this.numValues);

        double res = this.cumulValues.get(e).getValue();
        res -= this.cumulValues.get(b - 1).getValue();
        return res;
    }

    /**
     * Gets the Squared Mean for elements within index b..e, multiplied by weight. Note that
     * n*mean^2 = sum^2/n when mean := sum/n
     *
     * @param b index of begin element
     * @param e index of end element
     * @return the sum of squared mean
     */
    private double getSSM(int b, int e) {
        double res = this.getSumOfWeightedValues(b, e);
        return res * res / this.getSumOfWeights(b, e);
    }

    /**
     * Finds CB[i+completedRows] given that the result is at least
     * bp+(completedRows-1) and less than ep+(completedRows-1)
     * Complexity: O(ep-bp) <= O(m) @
     *
     * @param i startIndex
     * @param bp endindex
     * @param ep
     *
     * @return the index
     */
    private int findMaxBreakIndex(int i, int bp, int ep) {
        assert (bp < ep);
        assert (bp <= i);
        assert (ep <= i + 1);
        assert (i < this.bufferSize);
        assert (ep <= this.bufferSize);

        double minSSM = this.previousSSM[bp] + this.getSSM(bp + this.completedRows, i + this.completedRows);
        int foundP = bp;
        while (++bp < ep) {
            double currSSM = this.previousSSM[bp] + this.getSSM(bp + this.completedRows, i + this.completedRows);
            if (currSSM > minSSM) {
                minSSM = currSSM;

                foundP = bp;
            }
        }
        this.currentSSM[i] = minSSM;
        return foundP;
    }

    /**
     * Find CB[i+completedRows] for all i>=bi and i<ei given that the
     * results are at least bp+(completedRows-1) and less than
     * ep+(completedRows-1)
     * Complexity: O(log(ei-bi)*Max((ei-bi),(ep-bp)))
     * <= O(m*log(m))
     *
     *
     * @param bi
     * @param ei
     * @param bp
     * @param ep
     * @return
     */
    private void calcRange(int bi, int ei, int bp, int ep) {
        assert (bi <= ei);

        assert (ep <= ei);
        assert (bp <= bi);

        if (bi == ei) {
            return;
        }
        assert (bp < ep);

        int mi = (int) Math.floor((bi + ei) / 2);
        int mp = this.findMaxBreakIndex(mi, bp, Math.min(ep, mi + 1));

        assert (bp <= mp);
        assert (mp < ep);
        assert (mp <= mi);

        // solve first half of the sub-problems with lower 'half' of possible outcomes
        this.calcRange(bi, mi, bp, Math.min(mi, mp + 1));

        this.classBreaks[this.classBreaksIndex + mi] = mp; // store result for the middle element.

        // solve second half of the sub-problems with upper 'half' of possible outcomes
        this.calcRange(mi + 1, ei, mp, ep);
    }
    
    /**
     * Swaps the content of the two lists with each other.
     */
    private void swapArrays() {
        double[] temp = new double[previousSSM.length];
        System.arraycopy(previousSSM, 0, temp, 0, previousSSM.length);
        
        previousSSM = new double[currentSSM.length];
        System.arraycopy(currentSSM, 0, previousSSM, 0, currentSSM.length);
       
        currentSSM = new double[temp.length];
        System.arraycopy(temp, 0, currentSSM, 0, temp.length);
    }

    /**
     * Starting point of calculation of breaks.
     *
     * complexity: O(m*log(m)*k)
     */
    private void calcAll() {
        if (this.numBreaks >= 2) {
            this.classBreaksIndex = 0;
            for (this.completedRows = 1; this.completedRows < this.numBreaks - 1; ++this.completedRows) {
                this.calcRange(0, this.bufferSize, 0, this.bufferSize); // complexity: O(m*log(m))
                
                swapArrays();
                this.classBreaksIndex += this.bufferSize;
            }
        }
    }
    
    /**
     *  Does the internal processing to actually create the breaks.
     *
     * @param k number of breaks
     * @param vcpc asc ordered input of values and their occurence counts.
     */
    public static double[] classifyJenksFisherFromValueCountPairs(int k, List<ValueCountPair> vcpc){
        double[] breaksArray = new double[k];
        int m = vcpc.size();
        assert(k <= m); // PRECONDITION
        if (k==0)
            return breaksArray;
        JenksFisher jf = new JenksFisher(vcpc, k);
        if (k > 1) {
            // runs the actual calculation
            jf.calcAll();
            int lastClassBreakIndex = jf.findMaxBreakIndex(jf.bufferSize - 1, 0, jf.bufferSize);
            while (--k!=0) {
                // assign the break values to the result
                breaksArray[k]= vcpc.get(lastClassBreakIndex + k).getValue();
                assert(lastClassBreakIndex < jf.bufferSize);
                if (k > 1)
                {
                    jf.classBreaksIndex -= jf.bufferSize;
                    lastClassBreakIndex = jf.classBreaks[jf.classBreaksIndex + lastClassBreakIndex];
                }
            }
            assert(jf.classBreaks[jf.classBreaksIndex] ==jf.classBreaks[0]);
        }
        assert(k == 0);
        breaksArray[0] = vcpc.get(0).getValue(); // break for the first class is the minimum of the dataset.
        return breaksArray;
    }

    /**
     * Main entry point for creation of Jenks-Fisher natural breaks.
     *
     * @param values array of the values, do not need to be sorted.
     * @param k number of breaks to create
     * @return Array with breaks
     */
    public static List<Double> createJenksFisherBreaksArray(List<Double> values, int k){
        List<ValueCountPair>sortedUniqueValueCounts = getValueCountPairs(values);
        
        double[] breaksArray = null;
        if (sortedUniqueValueCounts.size()>k){
            breaksArray = classifyJenksFisherFromValueCountPairs(k, sortedUniqueValueCounts);
        }   
        else {
           breaksArray = new double[sortedUniqueValueCounts.size()];
           int i=0;
            for (ValueCountPair vcp : sortedUniqueValueCounts) {
                breaksArray[i] = vcp.getValue();
                i++;
            }

        }
        List<Double> result = new ArrayList<>(breaksArray.length);
        for(double d:breaksArray) result.add(d);
        return result;
    }

    /**
     * Calculates the occurence count of given values and returns them.
     *
     * @param values
     * @return Occurences of values.
     */
    private static List<ValueCountPair> getValueCountPairs(List<Double> values) {
        List<ValueCountPair> result = new ArrayList<>();
        HashMap<Double,ValueCountPair> vcpMap = new HashMap<>();
        for (double value : values) {
            if (!vcpMap.containsKey(value)) {
                ValueCountPair vcp = new ValueCountPair(value, 1);
                vcpMap.put(value, vcp);
                result.add(vcp);
            } else {
                vcpMap.get(value).incCount();
            }   
        }
        Collections.sort(result, new Comparator<ValueCountPair>(){
            @Override
            public int compare(ValueCountPair o1, ValueCountPair o2) {
                return Double.compare(o1.getValue(),o2.getValue());
            }
            
        });
        return result;
    }

    
}

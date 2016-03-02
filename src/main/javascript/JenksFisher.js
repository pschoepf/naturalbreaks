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

/**
 * Internal helper used for assertions.
 *
 * @param  {boolean} correct The boolean input.
 * @param  {String} msg     Message to be thrown.
 */
function assert(correct, msg) {
    //alert("correct="+correct+",\n\nactual="+actual);
    if (!correct)
        throw new Error("ASSERTION FAILED: " + msg);
}

/**
 * Internal helper to swap two arrays.
 *
 * @param  {Array} a Array one.
 * @param  {Array} b Array two and the holder of the result.
 */
function swapArrays(a, b) {
    var temp = a.slice(0);
    a.length = 0;
    a.push.apply(a, b);
    b.length = 0;
    b.push.apply(b, temp);
}

/**
* Constructor that initializes main variables used in fisher calculation of natural breaks.
*
* @param vcpc Ordered list of pairs of values to occurrence counts.
* @param k Number of breaks to find.
*/
function JenksFisher(/*array of {val:double, count:double}*/ vcpc, /*int*/ k) {
  if (!(this instanceof JenksFisher)) {
    return new JenksFisher();
  }
  this.cumulValues = [/*{wv: weighted value as double,cwv: cummul weighted value as double}*/];
  this.numValues = 0;
  this.numBreaks = 0;
  this.bufferSize = 0;
  this.previousSSM = [];
  this.currentSSM = [];
  this.classBreaks = [];
  this.classBreaksIndex = 0;
  this.completedRows = 0;
  this.cumulValues.length = [],
          this.numValues = vcpc.length;
  this.numBreaks = k,
          this.bufferSize = (vcpc.length - (k - 1)),
          this.previousSSM = new Array(this.bufferSize),
          this.currentSSM = new Array(this.bufferSize),
          this.classBreaks = new Array(this.bufferSize * (this.numBreaks - 1)),
          this.classBreaksIndex = 0, this.completedRows = 0;
  var cwv = 0.0;
  var cw = 0, w = 0;
  //console.debug("init:: values and weight:"+vcpc.length+", buffSize:"+this.bufferSize+", content:"+JSON.stringify(vcpc));
  for (var i = 0; i != this.numValues; ++i)
  {
      assert(!i || vcpc[i].val >= vcpc[i - 1].val); // PRECONDITION: the value sequence must be strictly increasing

      w = vcpc[i].count;
      assert(w > 0); // PRECONDITION: all weights must be positive

      cw += w;
      assert(cw >= w); // No overflow? No loss of precision?

      cwv += w * vcpc[i].val;
      this.cumulValues.push({cwv: cwv, cw: cw});
      if (i < this.bufferSize)
          this.previousSSM[i] = cwv * cwv / cw; // prepare SSM for first class. Last (k-1) values are omitted
  }
  //console.debug("init:: CumulValues:"+this.cumulValues);
  //console.debug("init:: previousSSM:"+this.previousSSM)
}


/**
 * Gets sum of weighs for elements with index b..e.
 *
 * @param b index of begin element
 * @param e index of end element
 * @return sum of weights.
 */
JenksFisher.prototype.getSumOfWeights = function (/*int*/ b, /*int*/ e)
        {
            assert(b);    // First element always belongs to class 0, thus queries should never include it.
            assert(b <= e);
            assert(e < this.numValues);

            var res = this.cumulValues[e].cw;
            res -= this.cumulValues[b - 1].cw;
            return res;
        }
/**
 * Gets sum of weighed values for elements with index b..e
 *
 * @param b index of begin element.
 * @param e index of end element
 * @return the cumul. sum of the values*weight
 */
JenksFisher.prototype.getSumOfWeightedValues = function (/*int*/ b, /*int*/ e)
        // Gets sum of weighed values for elements b..e
        {
            assert(b);
            assert(b <= e);
            assert(e < this.numValues);

            var res = this.cumulValues[e].cwv;
            res -= this.cumulValues[b - 1].cwv;
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
JenksFisher.prototype.getSSM = function (/*int*/ b, /*int*/ e)
{
    var res = this.getSumOfWeightedValues(b, e);
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
JenksFisher.prototype.findMaxBreakIndex = function (/*int*/ i, /*int*/ bp, /*int*/ ep)
{
    assert(bp < ep);
    assert(bp <= i);
    assert(ep <= i + 1);
    assert(i < this.bufferSize);
    assert(ep <= this.bufferSize);

    var minSSM = this.previousSSM[bp] + this.getSSM(bp + this.completedRows, i + this.completedRows);
    var foundP = bp;
    while (++bp < ep)
    {
        var currSSM = this.previousSSM[bp] + this.getSSM(bp + this.completedRows, i + this.completedRows);
        if (currSSM > minSSM)
        {
            //console.debug("findMaxBreakIndex:: new minSSM: " + minSSM+ ", currSSM:"+ currSSM);
            minSSM = currSSM;

            foundP = bp;
        }
    }
    this.currentSSM[i] = minSSM;
    return foundP;
},

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
JenksFisher.prototype.calcRange = function (/*int*/ bi, /*int*/ ei, /*int*/ bp, /*int*/ ep)
{
    assert(bi <= ei);

    assert(ep <= ei);
    assert(bp <= bi);

    if (bi === ei)
        return;
    assert(bp < ep);

    var mi = Math.floor((bi + ei) / 2);
    var mp = this.findMaxBreakIndex(mi, bp, Math.min(ep, mi + 1));

    assert(bp <= mp);
    assert(mp < ep);
    assert(mp <= mi);

    // solve first half of the sub-problems with lower 'half' of possible outcomes
    this.calcRange(bi, mi, bp, Math.min(mi, mp + 1));

    this.classBreaks[this.classBreaksIndex + mi ] = mp; // store result for the middle element.

    // solve second half of the sub-problems with upper 'half' of possible outcomes
    this.calcRange(mi + 1, ei, mp, ep);
},
/**
* Starting point of calculation of breaks.
*
* complexity: O(m*log(m)*k)
*/
JenksFisher.prototype.calcAll = function ()
    {
        if (this.numBreaks >= 2)
        {
            this.classBreaksIndex = 0;
            for (this.completedRows = 1; this.completedRows < this.numBreaks - 1; ++this.completedRows)
            {
                this.calcRange(0, this.bufferSize, 0, this.bufferSize); // complexity: O(m*log(m))

                swapArrays(this.previousSSM, this.currentSSM);
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
function classifyJenksFisherFromValueCountPairs(/*double[]*/ breaksArray, /*int*/ k, /*{val:value,count:count}*/ vcpc)
{
    var m = vcpc.length;
    assert(k <= m); // PRECONDITION
    if (!k)
        return;
    var jf = new JenksFisher(vcpc, k);
    if (k > 1) {
        // runs the actual calculation
        jf.calcAll();
        var lastClassBreakIndex = jf.findMaxBreakIndex(jf.bufferSize - 1, 0, jf.bufferSize);
        while (--k) {
            // assign the break values to the result
            breaksArray[k] = vcpc[lastClassBreakIndex + k].val;
            assert(lastClassBreakIndex < jf.bufferSize);
            if (k > 1)
            {
                jf.classBreaksIndex -= jf.bufferSize;
                lastClassBreakIndex = jf.classBreaks[jf.classBreaksIndex + lastClassBreakIndex];
            }
        }
        assert(jf.classBreaks[jf.classBreaksIndex] === jf.classBreaks[0]);
    }
    assert(k === 0);
    breaksArray[0] = vcpc[0].val; // break for the first class is the minimum of the dataset.
}

/**
 * Main entry point for creation of Jenks-Fisher natural breaks.
 *
 * @param {type} values - array of the values, do not need to be sorted.
 * @param {type} k - number of breaks to create
 * @returns {Array} - Array with breaks
 */
function createJenksFisherBreaksArray(/*double*/ values,/*int*/ k)
{
    // create pair of value->number of occurences (weight) which is input for the JF- algorithm
    var sortedUniqueValueCounts = getValueCountPairs(values);

    var breaksArray = new Array(k);
    if (sortedUniqueValueCounts.length>k){
         classifyJenksFisherFromValueCountPairs(breaksArray, k, sortedUniqueValueCounts);
    }
    else {
        breaksArray = [];
        sortedUniqueValueCounts.forEach(function(el){
            breaksArray.push(Math.abs(el.val))});
    }
    return breaksArray;
}
/**
 * create list with weighted unique values, sorted ASC by values.
 * @param {type} rawValues
 * @returns [{val:value as int,count:weight as int},...]
 */
function getValueCountPairs(rawValues) {
    var result = {};
    rawValues.forEach(function (el) {
        var k = "k_"+Number(el).toString();
        if (result[k]=="undefined"|| result[k] == null){
            result[k] = 1;
        }else{
            result[k] = result[k] + 1;
        }
    });
    var sortedResults = [];
    Object.keys(result).forEach(function (key) {
        sortedResults.push({val: Number(key.substring(2)), count: result[key]});
    });
    sortedResults.sort(function (a, b) {
        return a.val - b.val;
    });
    return sortedResults;
}

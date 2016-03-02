describe("Fisherbreaks javascript tests.", function() {
  it("should return expected fisher breaks results", function() {

    expect(createJenksFisherBreaksArray!=null).toBe(true);
    var input = [0.0,0.0,10.0,10.0,10.0,10.0,100.0,100.0,1002.0,10040.0,1009.0,101.0,101.0,1010.0,1011.0,1012.0,1015.0,1017.0,1018.0,102.0,102.0,1020.0];
    var expectedValues = [0.0,10.0,100.0,1002.0,1009.0,1011.0,1015.0,1017.0,1020.0,10040.0];
    var results = createJenksFisherBreaksArray(input, 10);
    expect(results).not.toEqual([0.47,0.74565,0.89474]);
    expect(results).toEqual(expectedValues);

  });

  it("should not take too long", function() {
      expect(createJenksFisherBreaksArray!=null).toBe(true);
      var input = [];
      for(var i=0;i<50000;i++) input.push(Math.random()*1000);
      var start = new Date().getTime();

      var results = createJenksFisherBreaksArray(input, 20);
      var timeTaken = (new Date().getTime()-start);

      expect(timeTaken<8000).toBe(true);
      it("performance test for n=50000,k=20 took in millis:" + timeTaken, function() {

         });
    });


});

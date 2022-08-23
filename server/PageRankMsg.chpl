module PageRankMsg {


  use Reflection;
  use ServerErrors;
  use Logging;
  use Message;
  use SegmentedArray;
  use ServerErrorStrings;
  use ServerConfig;
  use MultiTypeSymbolTable;
  use MultiTypeSymEntry;
  use RandArray;
  use IO;


  use SymArrayDmap;
  use Random;
  use RadixSortLSD;
  use Set;
  use DistributedBag;
  use ArgSortMsg;
  use Time;
  use CommAggregation;
  use Sort;
  use Map;
  use DistributedDeque;


  use List; 
  use LockFreeStack;
  use Atomics;
  use IO.FormattedIO; 
  use GraphArray;
  use GraphMsg;


  private config const logLevel = ServerConfig.logLevel;
  const smLogger = new Logger(logLevel);
  

  //Given a graph, calculate the pagerank of the graph
  proc segPageRankMsg(cmd: string, payload: string, st: borrowed SymTab): MsgTuple throws {
      var repMsg: string;
      var (n_verticesN,n_edgesN,directedN,weightedN,graphEntryName,restpart )
          = payload.splitMsgToTuple(6);
      var Nv=n_verticesN:int;
      var Ne=n_edgesN:int;
      var Directed=false:bool;
      var Weighted=false:bool;
      if (directedN:int)==1 {
          Directed=true;
      }
      if (weightedN:int)==1 {
          Weighted=true;
      }
      
      var countName:string;
      var timer:Timer;
      timer.start();

      var PageRank:[0..Nv-1] real;
      var tmpPageRank: [0..Nv-1] real;
      var subSum: [0..numLocales-1] real;
      var subError: [0..numLocales-1] real;
      var TotalSum: [0..0] real;
      var d = 0.85: real;
      var maxerror = 1.0 :real;
          
      PageRank=1.0/Nv;
      tmpPageRank=0.0;
      TotalSum = 0.0;
      subSum = 0.0;
      subError = 0.0;



      var gEntry:borrowed GraphSymEntry = getGraphSymEntry(graphEntryName, st);
      var ag = gEntry.graph;

      
      proc pageRank_kernel(nei:[?D1] int, start_i:[?D2] int,src:[?D3] int, dst:[?D4] int,
                        neiR:[?D11] int, start_iR:[?D12] int,srcR:[?D13] int, dstR:[?D14] int):string throws{

          
          var iteration = 1 :int;

          while(maxerror > 0.00001){
                maxerror = -1.0 : real;
                subSum=0.0;
                TotalSum=0.0;
                writeln(iteration, "round of iteration");
                iteration+=1;
                coforall loc in Locales {
                    on loc {
                        var ld = nei.localSubdomain();
                        var startVer = ld.low;
                        var endVer = ld.high;
                        var sum=0.0: real;
                  
                        forall i in startVer..endVer with (+ reduce sum){
                            if (neiR[i] != 0){
                                var beginTmp = start_iR[i];
                                var endTmp = beginTmp+neiR[i]-1;
                                forall j in beginTmp..endTmp {
                                    tmpPageRank[i]+=PageRank[dstR[j]]/nei[dstR[j]];
                                }
                            }

                            if (nei[i] == 0){
                                sum+=(PageRank[i]/(Nv:real));
                            }
                        }
                        subSum[here.id]=sum;
                    }
                    
                }

                for i in subSum{
                    TotalSum[0]+=i;
                }

                //writeln("TotalSum is ", TotalSum[0]);
                
                // for i in 0..Nv-1{
                //     tmpPageRank[i] = 0.85 * (tmpPageRank[i]+TotalSum[0]) + (1.0-0.85) / (Nv:real); 
                //     //writeln(i,"'s pagerank is ", tmpPageRank[i]);
                //     if (abs(PageRank[i]-tmpPageRank[i]) >= maxerror ) {
                //         maxerror = abs(PageRank[i]-tmpPageRank[i]);
                //     }
                // }

                coforall loc in Locales{
                    on loc{
                        var ld = nei.localSubdomain();
                        var startVer = ld.low;
                        var endVer = ld.high;
                        var subMaxError = -1.0:real;
                        forall i in startVer..endVer with (max reduce subMaxError){
                            tmpPageRank[i] = 0.85 * (tmpPageRank[i]+TotalSum[0]) + (1.0-0.85) / Nv:real;
                            //writeln(i,"'s pagerank is ", tmpPageRank[i]);
                            if (abs(PageRank[i]-tmpPageRank[i]) >= subMaxError ) {
                                subMaxError = abs(PageRank[i]-tmpPageRank[i]);
                            }
                        }
                        subError[here.id] = subMaxError;
                    }
                }

                for i in subError{
                    if (i > maxerror){
                        maxerror=i;
                    }
                }
                
                PageRank<=>tmpPageRank;
                tmpPageRank=0.0;

          }

          var countName = st.nextName();
          var countEntry = new shared SymEntry(PageRank);
          st.addEntry(countName, countEntry);

          var cntMsg =  'created ' + st.attrib(countName);
          return cntMsg;

      }





      if (!Directed) {
              repMsg=pageRank_kernel(
                      toSymEntry(ag.getNEIGHBOR(), int).a,
                      toSymEntry(ag.getSTART_IDX(), int).a,
                      toSymEntry(ag.getSRC(), int).a,
                      toSymEntry(ag.getDST(), int).a,
                      toSymEntry(ag.getNEIGHBOR_R(), int).a,
                      toSymEntry(ag.getSTART_IDX_R(), int).a,
                      toSymEntry(ag.getSRC_R(), int).a,
                      toSymEntry(ag.getDST_R(), int).a);
      }
      timer.stop();
      smLogger.debug(getModuleName(),getRoutineName(),getLineNumber(),repMsg);
      return new MsgTuple(repMsg, MsgType.NORMAL);
  }// end of seg







    proc registerMe() {
        use CommandMap;
        registerFunction("segmentedGraphPR", segPageRankMsg);
    }


}


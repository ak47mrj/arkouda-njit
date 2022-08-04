#!/usr/bin/env python3                                                         

import time, argparse
import numpy as np
import arkouda as ak
import random
import string
import arkouda_njit as njit

TYPES = ('int64', 'float64', 'bool', 'str')

def time_ak_test():
    print("Graph Truss Analysis")
    cfg = ak.get_config()
    print("server Hostname =",cfg["serverHostname"])
    print("Number of Locales=",cfg["numLocales"])
    print("number of PUs =",cfg["numPUs"])
    print("Max Tasks =",cfg["maxTaskPar"])
    print("Memory =",cfg["physicalMemory"])
    HomeDir="/rhome/zhihui/"
    Test1=[ \
            [3056,1024,2,0,HomeDir+"Adata/Delaunay/delaunay_n10/delaunay_n10.mtx.pr"],\
            [28980,5242,2,0,HomeDir+"Adata/SNAP/ca-GrQc.txt.gr.pr"],\
            [11,6,2,0,"g.gr"],\
              ]
    TestMtx=[ \

            [393176,131072,2,0,HomeDir+"Adata/Delaunay/delaunay_n17/delaunay_n17.mtx"],\
            [786396,262144,2,0,HomeDir+"Adata/Delaunay/delaunay_n18/delaunay_n18.mtx"],\
            [1572823,524288,2,0,HomeDir+"Adata/Delaunay/delaunay_n19/delaunay_n19.mtx"],\
            [3145686,1048576,2,0,HomeDir+"Adata/Delaunay/delaunay_n20/delaunay_n20.mtx"],\
            [6291408,2097152,2,0,HomeDir+"Adata/Delaunay/delaunay_n21/delaunay_n21.mtx"],\
            [12582869,4194304,2,0,HomeDir+"Adata/Delaunay/delaunay_n22/delaunay_n22.mtx"],\
            [25165784,8388608,2,0,HomeDir+"Adata/Delaunay/delaunay_n23/delaunay_n23.mtx"],\
            [50331601,16777216,2,0,HomeDir+"Adata/Delaunay/delaunay_n24/delaunay_n24.mtx"],\
              ]
    TestRGG=[ \
            [63501393,8388608,2,0,HomeDir+"Adata/rgg_n_2/rgg_n_2_23_s0/rgg_n_2_23_s0.mtx"],\
            [132557200,16777216,2,0,HomeDir+"Adata/rgg_n_2/rgg_n_2_24_s0/rgg_n_2_24_s0.mtx"],\
              ]
    TestKron=[ \
            [2456398,65536,3,0,HomeDir+"Adata/kron_g500-logn/kron_g500-logn16/kron_g500-logn16.mtx"],\
            [5114375,131072,3,0,HomeDir+"Adata/kron_g500-logn/kron_g500-logn17/kron_g500-logn17.mtx"],\
            [10583222,262144,3,0,HomeDir+"Adata/kron_g500-logn/kron_g500-logn18/kron_g500-logn18.mtx"],\
            [21781478,524288,3,0,HomeDir+"Adata/kron_g500-logn/kron_g500-logn19/kron_g500-logn19.mtx"],\
            [44620272,1048576,3,0,HomeDir+"Adata/kron_g500-logn/kron_g500-logn20/kron_g500-logn20.mtx"],\
              ]



    start = time.time()
    for i in Test1:
        Edges=i[0]
        Vertices=i[1]
        Columns=i[2]
        Directed=i[3]
        FileName=i[4]
        print(i)
        print(Edges,",",Vertices,",",Columns,",",Directed,",",str(FileName))
        Graph=njit.graph_file_read(Edges,Vertices,Columns,Directed,str(FileName),0,0,0,0,1)
        Jacc=njit.graph_jaccard_coefficient(Graph)
    end = time.time()
    return


def create_parser():
    parser = argparse.ArgumentParser(description="Measure the performance of suffix array building: C= suffix_array(V)")
    parser.add_argument('hostname', help='Hostname of arkouda server')
    parser.add_argument('port', type=int, help='Port of arkouda server')
    parser.add_argument('-v', '--logvertices', type=int, default=5, help='Problem size: log number of vertices')
    parser.add_argument('-e', '--vedges', type=int, default=2,help='Number of edges per vertex')
    parser.add_argument('-p', '--possibility', type=float, default=0.01,help='Possibility ')
    parser.add_argument('-t', '--trials', type=int, default=6, help='Number of times to run the benchmark')
    parser.add_argument('-m', '--perm', type=int, default=0 , help='if permutation ')
    parser.add_argument('--numpy', default=False, action='store_true', help='Run the same operation in NumPy to compare performance.')
    parser.add_argument('--correctness-only', default=False, action='store_true', help='Only check correctness, not performance.')
    return parser


    
if __name__ == "__main__":
    import sys
    
    parser = create_parser()
    args = parser.parse_args()
    ak.verbose = False
    print("Jaccard Coefficient Test")
    ak.connect(args.hostname, args.port)
    time_ak_test()

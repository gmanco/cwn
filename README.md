# CWN: An Expectation-Maximization framework for Influence-Based Network Oblivious Community Detection

CWN is a stochastic framework for detecting social communities when the social graph is not available, but instead we have access to a log of user activity, that is a dataset of tuples (u, i, t) recording the fact that user u “adopted” item i at time t.

The framework assumes that the adoption of items is governed by an underlying diffusion process over the unobserved social network, and that such diffusion model is based on community-level influence. By fitting the model parameters to the user activity log, we learn the community membership and the level of influence of each user in each community. 

The general framework is instantiated with two different diffusion models, one with discrete time and one with continuous time, and we show that the computational complexity of both approaches is linear in the number of users and in the size of the propagation log.

More details on the approach can be found at 

Nicola Barbieri, Francesco Bonchi, Giuseppe Manco. 2016. Influence-based Network-oblivious Community Detection. ACM Trans. Intell. Syst. Technol. 

- The package includes a synthetic dataset 
   - /resources/datasets/synthetic/s1/NR/fold1/actionLog

- The package 
    - src/propagationGenerator 

    contains classes to generate synthetic propagation given a network file

- The package 
    - src/topicBased 
    
    does clustering on the activation log (user that adopted the same items will be in the same cluster)

- The package
    - src/cicm
   
   instantiates CWN to the Community-based Independent Cascade Model CIC

- The package
    - src/diffusionBased
   
   instantiates CWN to the survival model CRate  

Each package contains a main class which can be executed. For example, If you want to run C-Rate:
- the main class is src/diffusionBased/CommunityRate_Inference.java

  taking parameters
  
   -a actionlog -c confFile -k nCommunities -o output -maxIt maxIt -g groundTruthCommunities  -l networkFile
  
  where:
   - actionLog is the file with node activations  (e.g. resources/datasets/synthetic/s1/NR/fold1/actionLog)
   - confFile describe the format of the actionLog (e.g. resources/datasets/synthetic/conf.inf)
   - maxIt specifies the maximum number of iterations
   - groundTruthCommunities is a mapping from node -> community (e.g. resources/datasets/synthetic/s1/NR/fold1/community.dat)
   - networkFile specifies the network (u \t v) => u follows v (e.g. resources/datasets/synthetic/s1/NR/fold1/network.dat)


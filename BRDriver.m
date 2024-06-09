function [ coding_result , lncRNA_result] = BRDriver( coding_tumor, lncRNA_tumor_profile)
%%   process data
tumor=importdata(coding_tumor);
tumor_data=tumor.data;

lncRNA_tumor=importdata(lncRNA_tumor_profile);
lncRNA_tumor_data=lncRNA_tumor.data;

patient_name = tumor.textdata(1,2:end);
coding_gene=tumor.textdata(2:end,1);
lncRNA_gene=lncRNA_tumor.textdata(2:end,1);
gene_list=[coding_gene;lncRNA_gene];
data=[tumor_data;lncRNA_tumor_data];            

%% Input PPI

Net_interactions = readtable('StringNet_lncRNA_mRNA.txt');
edge = table2array(Net_interactions);

[y1,~]=ismember(edge(:,1),gene_list);
[y2,~]=ismember(edge(:,2),gene_list);
y=y1.*y2;
new_edge = edge;
new_edge(y==0,:)=[];   
node = unique(new_edge);

[~,loc3]=ismember(node,gene_list);
new_gene_list = gene_list(loc3,:);
data = data(loc3,:);

[~,loc1]=ismember(new_edge(:,1),new_gene_list);
[~,loc2]=ismember(new_edge(:,2),new_gene_list);
location=[loc1 loc2];
[N2,~]=size(location);
N1=length(new_gene_list);
Net=zeros(N1,N1);
for i=1:N2
    Net(location(i,2),location(i,1))=1;
    Net(location(i,1),location(i,2))=1;
end

%% lioness  & RWR

coding_result = {};
lncRNA_result = {};
for i=1:size(data,2)
    
    tic
    i
    
    PN=lioness_method(data,i);
    
    adjacency_matrix=PN.*Net;
    
    %% Betweenness
    n = size(adjacency_matrix,1);
    Betweenness = betweenness_centrality( sparse(adjacency_matrix) );
    Betweenness = Betweenness/((n-1)*(n-2));
    Betweenness (Betweenness==0) = eps;
    
    %% RWR
    restart_probability = 0.85;
    node_importance = Random_walk_with_restart( adjacency_matrix,restart_probability );
    
    %% 
    MaxBetweenness = max( Betweenness(:,1) );
    Maxnode_importance = max( node_importance(:,1) );
    B = Betweenness/MaxBetweenness;
    RWR = node_importance/Maxnode_importance;
    Alpha = 0.5;
    Score = Alpha*B+(1-Alpha)*RWR;
    
    [ScoreRank,I] = sort( Score, 'descend');
    geneRank = new_gene_list(I);
    
    %%
    threshold = mean(Score); 
    num_nodes = size(Score,1);
    Index = 0;  
    for q=1:num_nodes
        if ScoreRank(q) > threshold
            Index = Index+1;
        else
            break
        end
    end
    
    Driver_gene = geneRank(1:Index,1);
    
    %%
    coding_driver = intersect(Driver_gene,coding_gene,'stable');
    lncRNA_driver = intersect(Driver_gene,lncRNA_gene,'stable');
    
    row = size(coding_driver,1);
    coding_result(2:row+1,i) = coding_driver ;  
    row = size(lncRNA_driver,1);
    lncRNA_result(2:row+1,i) = lncRNA_driver ;  
    
    toc
    
    clear  adjacency_matrix
end

coding_result(1,:) = patient_name;
lncRNA_result(1,:) = patient_name;

end
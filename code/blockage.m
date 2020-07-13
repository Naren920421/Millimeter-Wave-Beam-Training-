function blockCluster = blockage(T,N,blockDuration,blockProb)

% From t = 2 to t = T
blockCluster = zeros(N,T);
for t = 1:T-1 % except the first location
    if t == 1
        blockCluster(:,t) = binornd(1,blockProb,[N,1]); % random blockages to each cluster
        blockCluster(logical(blockCluster(:,t)==1),t+blockDuration-1) = 1;
    else
        blockCluster(logical(blockCluster(:,t)==0),t) = binornd(1,blockProb,[length(find(blockCluster(:,t)==0)),1]);
        if t > 1
            blockCluster(logical(blockCluster(:,t)==1&blockCluster(:,t-1)==0),t+blockDuration-1) = 1; % t>1, t<T
        end
        if t > 2
            blockCluster(logical(blockCluster(:,t-1)==1&blockCluster(:,t-2)==1),t) = 0; % t>2,t<=T
        end
        
    end
end



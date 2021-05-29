%%%  Visualize clusters for patterns of temporal gene expression
cd('F:\E Disk\3-科研合作\寄生虫课题\Hexa Network')
DEG=xlsread('DEG_AC Condition.xls');

for i=1:110
DEG(i,:)=DEG(i,:)/max(DEG(i,:));
end


%%% spine

DEG_spline=pchip([6 12 24],DEG,6:1:24);   %% 三次Hermite多项式插值，可用pchip函数替代

xlswrite('F:\寄生虫课题\Hexa Network\DEG_con.xls', DEG_spline);  

%%%%%%%%%  Clustering


% for i=1:110
%     for j=1:110
%         D(i,j)=pdist2(diff(DEG_spline(i,:)),diff(DEG_spline(j,:)),'hamming');
%     end
% end

% % Y=pdist(X,’metric’)      
% % 'euclidean'：欧氏距离（默认）；
% % 'seuclidean'：标准化欧氏距离；
% % 'mahalanobis'：马氏距离； 
% % 'cityblock'：布洛克距离；
% % 'minkowski'：明可夫斯基距离； 
% % 'cosine'： 夹角余弦
% % 'correlation'：  相关距离
% % 'spearman'             
% % 'hamming'： 汉明距离
% % 'jaccard'：    杰卡德距离& 杰卡德相似系数               
% % 'chebychev'：Chebychev距离

Y=pdist(DEG_spline,'cosine') ;    % 'cityblock'  
D = squareform(Y);
size(D)    
idx = kmedioids(D,6);

% idx = kmeans(DEG_spline,4);
% 
% idx = clusterdata(DEG_spline,5);



DEG_cluster_1_spline=DEG_spline(idx==1,:);
DEG_cluster_2_spline=DEG_spline(idx==2,:);
DEG_cluster_3_spline=DEG_spline(idx==3,:);
DEG_cluster_4_spline=DEG_spline(idx==4,:);
DEG_cluster_5_spline=DEG_spline(idx==5,:);
DEG_cluster_6_spline=DEG_spline(idx==6,:);
% DEG_cluster_7_spline=DEG_spline(idx==7,:);
% DEG_cluster_8_spline=DEG_spline(idx==8,:);
% sum(idx==6)
% size(DEG_cluster_1_spline)
 
 

close all

%%%% cluster 1
figure,
plot([6:1:24],DEG_cluster_1_spline,'Color',[0.5 0.5 0.5],'MarkerSize',2);axis([6 24 0 1.5]);
set(gca,'xtick',[6 12 24],'FontSize',16)
set(gca,'xticklabel',{'6','12','24'})
xlabel('Time (Days)','FontSize',18);
ylabel('Gene Expression','FontSize',18);
hold on,
plot([6:1:24],mean(DEG_cluster_1_spline,1),'r-','LineWidth',2);
hold on,
plot([6 12 24],mean(DEG_cluster_1_spline(:,[6,12,24]-5),1),'r.','MarkerSize',30);
title('cluster 1')

%%%% cluster 2
figure,
plot([6:1:24],DEG_cluster_2_spline,'Color',[0.5 0.5 0.5],'MarkerSize',2);axis([6 24 0 1.5]);
set(gca,'xtick',[6 12 24],'FontSize',16)
set(gca,'xticklabel',{'6','12','24'})
xlabel('Time (Days)','FontSize',18);
ylabel('Gene Expression','FontSize',18);
hold on,
plot([6:1:24],mean(DEG_cluster_2_spline,1),'r-','LineWidth',2);
hold on,
plot([6 12 24],mean(DEG_cluster_2_spline(:,[6,12,24]-5),1),'r.','MarkerSize',30);
title('cluster 2')


%%%% cluster 3
figure,
plot([6:1:24],DEG_cluster_3_spline,'Color',[0.5 0.5 0.5],'MarkerSize',2);axis([6 24 0 1.5]);
set(gca,'xtick',[6 12 24],'FontSize',16)
set(gca,'xticklabel',{'6','12','24'})
xlabel('Time (Days)','FontSize',18);
ylabel('Gene Expression','FontSize',18);
hold on,
plot([6:1:24],mean(DEG_cluster_3_spline,1),'r-','LineWidth',2);
hold on,
plot([6 12 24],mean(DEG_cluster_3_spline(:,[6,12,24]-5),1),'r.','MarkerSize',30);
title('cluster 3')


%%%% cluster 4
figure,
plot([6:1:24],DEG_cluster_4_spline,'Color',[0.5 0.5 0.5],'MarkerSize',2);axis([6 24 0 1.5]);
set(gca,'xtick',[6 12 24],'FontSize',16)
set(gca,'xticklabel',{'6','12','24'})
xlabel('Time (Days)','FontSize',18);
ylabel('Gene Expression','FontSize',18);
hold on,
plot([6:1:24],mean(DEG_cluster_4_spline,1),'r-','LineWidth',2);
hold on,
plot([6 12 24],mean(DEG_cluster_4_spline(:,[6,12,24]-5),1),'r.','MarkerSize',30);
title('cluster 4')


%%%% cluster 5
figure,
plot([6:1:24],DEG_cluster_5_spline','Color',[0.5 0.5 0.5],'MarkerSize',2);axis([6 24 0 1.5]);
set(gca,'xtick',[6 12 24],'FontSize',16)
set(gca,'xticklabel',{'6','12','24'})
xlabel('Time (Days)','FontSize',18);
ylabel('Gene Expression','FontSize',18);
hold on,
plot([6:1:24],mean(DEG_cluster_5_spline,1),'r-','LineWidth',2);
hold on,
plot([6 12 24],mean(DEG_cluster_5_spline(:,[6,12,24]-5),1),'r.','MarkerSize',30);
title('cluster 5')


%%%% cluster 6
figure,
plot([6:1:24],DEG_cluster_6_spline,'Color',[0.5 0.5 0.5],'MarkerSize',2);axis([6 24 0 1.5]);
set(gca,'xtick',[6 12 24],'FontSize',16)
set(gca,'xticklabel',{'6','12','24'})
xlabel('Time (Days)','FontSize',18);
ylabel('Gene Expression','FontSize',18);
hold on,
plot([6:1:24],mean(DEG_cluster_6_spline,1),'r-','LineWidth',2);
hold on,
plot([6 12 24],mean(DEG_cluster_6_spline(:,[6,12,24]-5),1),'r.','MarkerSize',30);
title('cluster 6')


% 
% 
% %%%% cluster 7
% figure,
% plot([6:1:24],DEG_cluster_7_spline','Color',[0.5 0.5 0.5],'MarkerSize',2);axis([6 24 0 1.5]);
% set(gca,'xtick',[6 12 24],'FontSize',16)
% set(gca,'xticklabel',{'6','12','24'})
% xlabel('Time (Days)','FontSize',18);
% ylabel('Gene Expression','FontSize',18);
% hold on,
% plot([6:1:24],mean(DEG_cluster_7_spline,1),'r-','LineWidth',2);
% hold on,
% plot([6 12 24],mean(DEG_cluster_7_spline(:,[6,12,24]-5),1),'r.','MarkerSize',30);
% title('cluster 7')
% 
% 
% %%%% cluster 8
% figure,
% plot([6:1:24],DEG_cluster_8_spline,'Color',[0.5 0.5 0.5],'MarkerSize',2);axis([6 24 0 1.5]);
% set(gca,'xtick',[6 12 24],'FontSize',16)
% set(gca,'xticklabel',{'6','12','24'})
% xlabel('Time (Days)','FontSize',18);
% ylabel('Gene Expression','FontSize',18);
% hold on,
% plot([6:1:24],mean(DEG_cluster_8_spline,1),'r-','LineWidth',2);
% hold on,
% plot([6 12 24],mean(DEG_cluster_8_spline(:,[6,12,24]-5),1),'r.','MarkerSize',30);
% title('cluster 8')
% 
% 

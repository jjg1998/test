clear
close all
clc
disp('The only input needed is a matrix file')
disp('The size of the data matrix is: N*d')
disp('N is the number of data points, d is the dimension of data points')

%%从文件中读取数据
mdist=input('name of the distance matrix file\n','s');
disp('Reading input matrix')
xx=load(mdist);%载入文件数据
[N,d]=size(xx);%N为数据点的总个数，d为数据点的维数

%输入所搜寻的k近邻的k值
disp("Please enter the needed k value")
k=input('k=');

%将输入矩阵转换为DPC算法所需矩阵形式
mdist=pdist(xx);         %两两行之间距离
A=tril(ones(N))-eye(N);
[x,y]=find(A~=0);
%第一列：元素i的标号，第二列：元素j的标号，第三列：元素i与元素j的距离
xxx=[x y mdist'];
ND=size(xxx,1); %% xx 第一个维度的长度，相当于文件的行数（即距离的总个数）

%计算距离矩阵dist(i,j)
%初始化为0
dist=zeros(N,N);
%进行赋值
for i=1:ND
  ii=xxx(i,1);
  jj=xxx(i,2);
  dist(ii,jj)=xxx(i,3);
  dist(jj,ii)=xxx(i,3);
end

%确定dc
percent=0.5;
fprintf('average percentage of neighbours (hard coded):%5.6f\n',percent);
position=round(ND*percent/100);%round是一个四舍五入函数
sda=sort(xxx(:,3));%对所有的距离值做升序排序
dc=sda(position);

%计算局部密度rho（利用Gaussian核）
fprintf('Computing Rho with gaussian kernel of radius:%12.6f\n',dc);
%将每个数据点的 rho 值初始化为零
rho=zeros(1,N);
%高斯核计算
for i=1:N-1
  for j=i+1:N
     rho(i)=rho(i)+exp(-(dist(i,j)/dc)*(dist(i,j)/dc));
     rho(j)=rho(j)+exp(-(dist(i,j)/dc)*(dist(i,j)/dc));
  end
end

%先求矩阵列最大值，再求最大值，最后得到所有距离值中的最大值
maxd=max(max(dist));
%将 rho 按降序排列，ordrho 保持序
[rho_sorted,ordrho]=sort(rho,'descend');
%处理 rho 值最大的数据点
delta(ordrho(1))=-1.;
nneigh(ordrho(1))=0;
%生成 delta 和 nneigh 数组
for ii=2:N
   delta(ordrho(ii))=maxd;
   for jj=1:ii-1
     if(dist(ordrho(ii),ordrho(jj))<delta(ordrho(ii)))
        delta(ordrho(ii))=dist(ordrho(ii),ordrho(jj));
        nneigh(ordrho(ii))=ordrho(jj);
        % 记录 rho 值更大的数据点中与 ordrho(ii) 距离最近的点的编号 ordrho(jj)
     end
   end
end
% 生成 rho 值最大数据点的 delta 值
delta(ordrho(1))=max(delta(:));

%取ind为下标索引，定义参数gamma=rho*delta,theta=delta/rho
%初始化为0
ind=zeros(1,N);
gamma=zeros(1,N);
theta=zeros(1,N);
%进行赋值
for i=1:N
  ind(i)=i;%下标
  gamma(i)=rho(i)*delta(i);
  theta(i)=delta(i)/rho(i);
end
%将gamma以及theta进行降序排序，ordgamma、ordtheta为保持序
[gamma_sorted,ordgamma]=sort(gamma,'descend');
[theta_sorted,ordtheta]=sort(theta,'descend');

%决策图1
disp('Generated file:DECISION GRAPH 1')
disp('column 1:ind')
disp('column 2:sorted gamma')
fid = fopen('DECISION_GRAPH_1', 'w');
for i=1:N
   fprintf(fid, '%6.2f %6.2f\n',i,gamma_sorted(i));
end

%每台计算机，句柄的根对象只有一个，就是屏幕，它的句柄总是 0
%>> scrsz = get(0,'ScreenSize')
%scrsz =
%           1           1        1280         800
%1280 和 800 就是你设置的计算机的分辨率，scrsz(4) 就是 800，scrsz(3) 就是 1280
scrsz = get(0,'ScreenSize');%ScreenSize is 四维向量: [left, bottom, width, height]

%人为指定一个位置
figure('Position',[6 72 scrsz(3)/2 scrsz(4)/1.3]);

%选择一个围住类中心的矩形
disp('Select a rectangle enclosing cluster centers')

%画出所谓的“决策图1”
subplot(1,1,1)
% tt=plot(1:N,gamma_sorted(:),'o','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k');
plot(rho(:),delta(:),'o','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k');
%o：圆圈；标记点大小：5；标记点内部填充颜色：黑色；标记点边框颜色：黑色；
title ('Decision Graph 1','FontSize',10.0)%画图中的text属性
xlabel ('IND')
ylabel ('Sorted gamma')

fig=subplot(1,1,1);
rect=getrect(fig);
%getrect 从图中用鼠标截取一个矩形区域， rect 中存放的是
%矩形左下角的坐标 (x,y) 以及所截矩形的宽度和高度
% gammamin=rect(2);
rhomin=rect(1);
deltamin=rect(2);

%找到每个数据点的k近邻
knn=zeros(N,k);
for i=1:N
    for j=1:k
        [dist_i_sorted,orddist_i]=sort(dist(i,:));
        knn(i,j)=orddist_i(j+1);%由于对角线元素均为0，所以要加1
    end
end

%定义点间连通性，不考虑对角线元素,对角线元素为0
connect=zeros(N,N);
for i=1:N-1
    for j=i+1:N
            if any(i==knn(j,:))
%any函数完成判断i是否等于knn的第j行的某一个元素。
                connect(i,j)=1;
                connect(j,i)=1;
                continue;
            end
            jn=j;
            jn_1=j;
            while connect(i,j)==0
                min=dist(i,jn_1);%在jn_1的k近邻内寻找与点i距离最近的点jn
                for ii=1:k
                    if min>dist(i,knn(jn_1,ii))
                        min=dist(i,knn(jn_1,ii));
                        jn=knn(jn_1,ii);
                    end
                end
                if any(i==knn(jn,:))
                    connect(i,j)=1;
                    connect(j,i)=1;
                    break;
                end
                if jn_1==jn
%点jn没有变化，表示点jn_1的k近邻内所有点与点i的距离，都大于等于点jn_1与点i的距离。
                    connect(i,j)=-1;
                    connect(j,i)=-1;
                    break;
                end
                jn_1=jn;
            end
    end
end

%初始化 cluster 个数
NCLUST=0;
%cl 为归属标志数组，cl(i)=j 表示第 i 号数据点归属于第 j 个 cluster
%先统一将 cl 初始化为 -1
cl=zeros(1,N);
for i=1:N
  cl(i)=-1;
end
for i=1:N
%   if gamma(i)>gammamin
   if ( (rho(i)>rhomin) && (delta(i)>deltamin))
     NCLUST=NCLUST+1;
     cl(i)=NCLUST;  %第 i 号数据点属于第 NCLUST 个 cluster
     icl(NCLUST)=i; %逆映射,第 NCLUST 个 cluster 的中心为第 i 号数据点
  end
end
fprintf('NUMBER OF CLUSTERS: %i \n', NCLUST);

if d==2
cmap=colormap;%返回当前图窗的颜色图，为RGB三元组成的3列矩阵,为渐变色。
disp('performing 2D scaling')
figure
plot(xx(:,1),xx(:,2),'o','MarkerSize',2,'MarkerFaceColor','k','MarkerEdgeColor','k');
hold on
for i=1:NCLUST
    ic=int8((i*64.)/(NCLUST*1.));
    for j=1:N
        if (cl(j)==i)
            plot(xx(j,1),xx(j,2),'o','MarkerSize',6,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
            hold on
        end
    end
end
title ('2D scaling','FontSize',10.0)
xlabel ('X')
ylabel ('Y') 
end

%新建空数组Aq与Bq，并把聚类中心依次放入Aq数组
Aq=[];
Bq=[];
for i=1:NCLUST
    Aq(i)=ordgamma(i);
end

%定义密度峰值点（若非聚类中心点的指向点的k近邻内，不包含该点，则称该点为密度峰值点）
%将所有与指向点连通，并且指向点也已经分配的密度峰值点加入Aq
%将与指向点不连通，或将与指向点连通但指向点并未分配的密度峰值点加入Bq
for i=1:N  
%按gamma降序排列对于每一个数据点操作
    if (cl(ordgamma(i))>0||any(ordgamma(i)==knn(nneigh(ordgamma(i)),:)))
%除去聚类中心点以及非密度峰值点(即点i在其指向点的k近邻内）
        continue;
    else
        if (connect(ordgamma(i),nneigh(ordgamma(i)))==1&&cl(nneigh(ordgamma(i)))>0)
%将与指向点连通，并且指向点也已分配的密度峰值点加入Aq
            cl(ordgamma(i))=cl(nneigh(ordgamma(i)));
            Aq=[Aq ordgamma(i)];
        else
%将与指向点不连通，或是与指向点连通但指向点并未分配的密度峰值点加入Bq
            Bq=[Bq ordgamma(i)];
        end
    end
end

%将Bq中所有与Aq存在连通的点依次加入Aq
while ~isempty(Bq)  %当Bq矩阵不为空时
    [i,j,Bq_j,max_d]=mindist(Aq,Bq,dist);
%该子函数找到Aq数组内点i与Bq数组内点j，使两点间距离最小，
%并输出Bq中点的所在位置Bq_j,以及Aq，Bq数组的最大距离
    if connect(i,j)==1  %当Aq与Bq中距离最近的两点是连通的
        cl(j)=cl(i);
        Aq=[Aq j];
        Bq(:,[Bq_j])=[];%删除Bq中的点
    else  %点i与点j不连通，将其距离改为最大距离
        dist(i,j)=max_d;
    end
    if isempty(Bq)
        break;
    end
    [logic]=link(Aq,Bq,connect);
    if logic==0
        %当Aq数组中的点与Bq数组内各点均不具有连通性时，退出循环
        break;
    end
end
        
%经上述处理后，若Bq不为空，此时Bq内各点均不与Aq内的点具有连通性
while ~isempty(Bq)
%找到Aq数组内点i与Bq数组内点j，使两点间距离最小，并输出Bq中点的所在位置Bq_j
    [i,j,Bq_j,~]=mindist(Aq,Bq,dist);
    cl(j)=cl(i);
    Bq(:,[Bq_j])=[];%删除Bq中的点
    Aq=[Aq j];%将该点加入Aq末尾
    if isempty(Bq)
        break;
    end
end

while any(-1==cl(:))    
for i=1:N
    if cl(ordgamma(i))==-1&&cl(nneigh(ordgamma(i)))>0;
        cl(ordgamma(i))=cl(nneigh(ordgamma(i)));
    end
end
end
%到此为止，完成聚类分配   

%如果为二维数据集，绘制聚类结果
if d==2
cmap=colormap;%返回当前图窗的颜色图，为RGB三元组成的3列矩阵,为渐变色。
disp('performing 2D clustering')
figure
for i=1:NCLUST
    ic=int8((i*64.)/(NCLUST*1.));
    for j=1:N
        if (cl(j)==i)
            plot(xx(j,1),xx(j,2),'o','MarkerSize',2,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
            hold on
        end
    end
end
    title ('2D Clustering','FontSize',10.0)
    xlabel ('X')
    ylabel ('Y')
end

%逐一处理每个 cluster
for i=1:NCLUST
  nc=0;  %% 用于累计当前 cluster 中数据点的个数
  for j=1:N
    if (cl(j)==i)
      nc=nc+1;
    end
  end
  fprintf('CLUSTER: %i CENTER: %i ELEMENTS: %i \n', i,icl(i),nc);
end

%存储分配结果
faa=fopen('CLUSTER_ASSIGNATION','w');
disp('Generated file:CLUSTER_ASSIGNATION')
disp('column 1:element id')
disp('column 2:cluster assignation without noise control')
for i=1:N
   fprintf(faa, '%i %i\n',i,cl(i));
end

%% 进行噪点识别
disp('Do you want to do noise recognition')
disp("'y' for yes,'n' for no.Please input")
answer=input('answer=');

if answer=='y'
disp('Generated file:DECISION GRAPH 2')
disp('column 1:IND')
disp('column 2:THETA SORTED')
 
fid = fopen('DECISION_GRAPH_2', 'w');
for i=1:N
   fprintf(fid, '%6.2f %6.2f\n', ind(i),theta_sorted(i));
end

%选择一个可以围住类噪点的矩形
disp('Select a rectangle enclosing cluster noises')

%画出所谓的“决策图2”
%人为指定一个位置
figure('Position',[6 72 scrsz(3)/2 scrsz(4)/1.3]);
subplot(1,1,1)
tt=plot(1:N,theta_sorted(:),'o','MarkerSize',3,'MarkerFaceColor','k','MarkerEdgeColor','k');
%o：圆圈；标记点大小：3；标记点内部填充颜色：黑色；标记点边框颜色：黑色；
title ('Decision Graph 2','FontSize',10.0)%画图中的text属性
xlabel ('IND')
ylabel ('Sorted \theta')

fig=subplot(1,1,1);
rect = getrect(fig);
%getrect 从图中用鼠标截取一个矩形区域， rect 中存放的是
%矩形左下角的坐标 (x,y) 以及所截矩形的宽度和高度

thetamin=rect(2);

%统计所划分区域的噪点个数
N_NOISE=0;
for i=1:N
  if ( theta(i)>thetamin)
     N_NOISE=N_NOISE+1;
  end
end

%令noise向量内部存储噪点的标号
for i=1:N_NOISE
    noise(i)=ordtheta(i);
end

%halo(i)=0代表该点为噪声点
for i=1:N
    halo(i)=cl(i);
    for j=1:N_NOISE
    if i==noise(j)
        halo(i)=0;
    end
    end
end
%调整noise向量中的顺序，以便于后面比较
noise=find(halo==0);

%参照k近邻，使用k近邻优化噪声识别
%初始化nrho(i)=0,定义噪点密度
nrho=zeros(1,N);
%nrho(i)>0表示点i属于核心点，nrho(i)<=0表示点i属于噪点
for i=1:N
    for j=1:k
            if halo(knn(i,j))>0
                nrho(i)=nrho(i)+1;
            else
                nrho(i)=nrho(i)-1;
            end
    end
end
%存储有可能成为噪点的点（取nrho<=k-1的数据点）
may_noise=find(nrho<=k-1);

while 1
N_NOISE=size(noise,2);    
for i=1:size(may_noise,2)
    for j=1:k
            if halo(knn(may_noise(i),j))>0
                nrho(may_noise(i))=nrho(may_noise(i))+1;
            else
                nrho(may_noise(i))=nrho(may_noise(i))-1;
            end
    end
end
%判断是否为噪点,nrho(i)<=0表示点i属于噪点
for i=1:N
    if nrho(i)<=0
        halo(i)=0;%点i为噪点
    else
        halo(i)=cl(i);%点i不为噪点
    end
end

NOISE=find(halo==0);

%判断NOISE向量是否发生变化，若未发生变化则说明噪点不在变化，退出循环
if N_NOISE==size(NOISE,2)
    if noise==NOISE
        break;
    end
end
noise=NOISE;
end
fprintf('NUMBER OF NOISES: %i \n', N_NOISE);

if d==2
%画出不带有噪声点的聚类结果
cmap=colormap;%返回当前图窗的颜色图，为RGB三元组成的3列矩阵,为渐变色。
disp('performing 2D clustering with noise control')
figure
for i=1:NCLUST
    ic=int8((i*64.)/(NCLUST*1.));
    for j=1:N
        if (halo(j)==i)
            plot(xx(j,1),xx(j,2),'o','MarkerSize',2,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
            hold on
        end
    end
end
title ('2D Clustering with noise control','FontSize',10.0)
xlabel ('X')
ylabel ('Y')
end

%逐一处理每个 cluster
for i=1:NCLUST
  nc=0;  %% 用于累计当前 cluster 中数据点的个数
  nh=0;  %% 用于累计当前 cluster 中核心数据点的个数
  for j=1:N
    if (cl(j)==i)
      nc=nc+1;
    end
    if (halo(j)==i)
      nh=nh+1;
    end
  end
  fprintf('CLUSTER: %i CENTER: %i ELEMENTS: %i CORE: %i NOISE: %i \n', i,icl(i),nc,nh,nc-nh);
end

%存储分配结果
faa=fopen('CLUSTER_ASSIGNATION_NOISE','w');
disp('Generated file:CLUSTER_ASSIGNATION_NOISE')
disp('column 1:element id')
disp('column 2:cluster assignation without noise control')
disp('column 3:cluster assignation with noise control')
for i=1:N
   fprintf(faa, '%i %i %i\n',i,cl(i),halo(i));
end
end
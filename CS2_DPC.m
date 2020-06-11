clear
close all
clc
disp('The only input needed is a matrix file')
disp('The size of the data matrix is: N*d')
disp('N is the number of data points, d is the dimension of data points')

%%���ļ��ж�ȡ����
mdist=input('name of the distance matrix file\n','s');
disp('Reading input matrix')
xx=load(mdist);%�����ļ�����
[N,d]=size(xx);%NΪ���ݵ���ܸ�����dΪ���ݵ��ά��

%��������Ѱ��k���ڵ�kֵ
disp("Please enter the needed k value")
k=input('k=');

%���������ת��ΪDPC�㷨���������ʽ
mdist=pdist(xx);         %������֮�����
A=tril(ones(N))-eye(N);
[x,y]=find(A~=0);
%��һ�У�Ԫ��i�ı�ţ��ڶ��У�Ԫ��j�ı�ţ������У�Ԫ��i��Ԫ��j�ľ���
xxx=[x y mdist'];
ND=size(xxx,1); %% xx ��һ��ά�ȵĳ��ȣ��൱���ļ�����������������ܸ�����

%����������dist(i,j)
%��ʼ��Ϊ0
dist=zeros(N,N);
%���и�ֵ
for i=1:ND
  ii=xxx(i,1);
  jj=xxx(i,2);
  dist(ii,jj)=xxx(i,3);
  dist(jj,ii)=xxx(i,3);
end

%ȷ��dc
percent=0.5;
fprintf('average percentage of neighbours (hard coded):%5.6f\n',percent);
position=round(ND*percent/100);%round��һ���������뺯��
sda=sort(xxx(:,3));%�����еľ���ֵ����������
dc=sda(position);

%����ֲ��ܶ�rho������Gaussian�ˣ�
fprintf('Computing Rho with gaussian kernel of radius:%12.6f\n',dc);
%��ÿ�����ݵ�� rho ֵ��ʼ��Ϊ��
rho=zeros(1,N);
%��˹�˼���
for i=1:N-1
  for j=i+1:N
     rho(i)=rho(i)+exp(-(dist(i,j)/dc)*(dist(i,j)/dc));
     rho(j)=rho(j)+exp(-(dist(i,j)/dc)*(dist(i,j)/dc));
  end
end

%������������ֵ���������ֵ�����õ����о���ֵ�е����ֵ
maxd=max(max(dist));
%�� rho ���������У�ordrho ������
[rho_sorted,ordrho]=sort(rho,'descend');
%���� rho ֵ�������ݵ�
delta(ordrho(1))=-1.;
nneigh(ordrho(1))=0;
%���� delta �� nneigh ����
for ii=2:N
   delta(ordrho(ii))=maxd;
   for jj=1:ii-1
     if(dist(ordrho(ii),ordrho(jj))<delta(ordrho(ii)))
        delta(ordrho(ii))=dist(ordrho(ii),ordrho(jj));
        nneigh(ordrho(ii))=ordrho(jj);
        % ��¼ rho ֵ��������ݵ����� ordrho(ii) ��������ĵ�ı�� ordrho(jj)
     end
   end
end
% ���� rho ֵ������ݵ�� delta ֵ
delta(ordrho(1))=max(delta(:));

%ȡindΪ�±��������������gamma=rho*delta,theta=delta/rho
%��ʼ��Ϊ0
ind=zeros(1,N);
gamma=zeros(1,N);
theta=zeros(1,N);
%���и�ֵ
for i=1:N
  ind(i)=i;%�±�
  gamma(i)=rho(i)*delta(i);
  theta(i)=delta(i)/rho(i);
end
%��gamma�Լ�theta���н�������ordgamma��ordthetaΪ������
[gamma_sorted,ordgamma]=sort(gamma,'descend');
[theta_sorted,ordtheta]=sort(theta,'descend');

%����ͼ1
disp('Generated file:DECISION GRAPH 1')
disp('column 1:ind')
disp('column 2:sorted gamma')
fid = fopen('DECISION_GRAPH_1', 'w');
for i=1:N
   fprintf(fid, '%6.2f %6.2f\n',i,gamma_sorted(i));
end

%ÿ̨�����������ĸ�����ֻ��һ����������Ļ�����ľ������ 0
%>> scrsz = get(0,'ScreenSize')
%scrsz =
%           1           1        1280         800
%1280 �� 800 ���������õļ�����ķֱ��ʣ�scrsz(4) ���� 800��scrsz(3) ���� 1280
scrsz = get(0,'ScreenSize');%ScreenSize is ��ά����: [left, bottom, width, height]

%��Ϊָ��һ��λ��
figure('Position',[6 72 scrsz(3)/2 scrsz(4)/1.3]);

%ѡ��һ��Χס�����ĵľ���
disp('Select a rectangle enclosing cluster centers')

%������ν�ġ�����ͼ1��
subplot(1,1,1)
% tt=plot(1:N,gamma_sorted(:),'o','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k');
plot(rho(:),delta(:),'o','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k');
%o��ԲȦ����ǵ��С��5����ǵ��ڲ������ɫ����ɫ����ǵ�߿���ɫ����ɫ��
title ('Decision Graph 1','FontSize',10.0)%��ͼ�е�text����
xlabel ('IND')
ylabel ('Sorted gamma')

fig=subplot(1,1,1);
rect=getrect(fig);
%getrect ��ͼ��������ȡһ���������� rect �д�ŵ���
%�������½ǵ����� (x,y) �Լ����ؾ��εĿ�Ⱥ͸߶�
% gammamin=rect(2);
rhomin=rect(1);
deltamin=rect(2);

%�ҵ�ÿ�����ݵ��k����
knn=zeros(N,k);
for i=1:N
    for j=1:k
        [dist_i_sorted,orddist_i]=sort(dist(i,:));
        knn(i,j)=orddist_i(j+1);%���ڶԽ���Ԫ�ؾ�Ϊ0������Ҫ��1
    end
end

%��������ͨ�ԣ������ǶԽ���Ԫ��,�Խ���Ԫ��Ϊ0
connect=zeros(N,N);
for i=1:N-1
    for j=i+1:N
            if any(i==knn(j,:))
%any��������ж�i�Ƿ����knn�ĵ�j�е�ĳһ��Ԫ�ء�
                connect(i,j)=1;
                connect(j,i)=1;
                continue;
            end
            jn=j;
            jn_1=j;
            while connect(i,j)==0
                min=dist(i,jn_1);%��jn_1��k������Ѱ�����i��������ĵ�jn
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
%��jnû�б仯����ʾ��jn_1��k���������е����i�ľ��룬�����ڵ��ڵ�jn_1���i�ľ��롣
                    connect(i,j)=-1;
                    connect(j,i)=-1;
                    break;
                end
                jn_1=jn;
            end
    end
end

%��ʼ�� cluster ����
NCLUST=0;
%cl Ϊ������־���飬cl(i)=j ��ʾ�� i �����ݵ�����ڵ� j �� cluster
%��ͳһ�� cl ��ʼ��Ϊ -1
cl=zeros(1,N);
for i=1:N
  cl(i)=-1;
end
for i=1:N
%   if gamma(i)>gammamin
   if ( (rho(i)>rhomin) && (delta(i)>deltamin))
     NCLUST=NCLUST+1;
     cl(i)=NCLUST;  %�� i �����ݵ����ڵ� NCLUST �� cluster
     icl(NCLUST)=i; %��ӳ��,�� NCLUST �� cluster ������Ϊ�� i �����ݵ�
  end
end
fprintf('NUMBER OF CLUSTERS: %i \n', NCLUST);

if d==2
cmap=colormap;%���ص�ǰͼ������ɫͼ��ΪRGB��Ԫ��ɵ�3�о���,Ϊ����ɫ��
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

%�½�������Aq��Bq�����Ѿ����������η���Aq����
Aq=[];
Bq=[];
for i=1:NCLUST
    Aq(i)=ordgamma(i);
end

%�����ܶȷ�ֵ�㣨���Ǿ������ĵ��ָ����k�����ڣ��������õ㣬��Ƹõ�Ϊ�ܶȷ�ֵ�㣩
%��������ָ�����ͨ������ָ���Ҳ�Ѿ�������ܶȷ�ֵ�����Aq
%����ָ��㲻��ͨ������ָ�����ͨ��ָ��㲢δ������ܶȷ�ֵ�����Bq
for i=1:N  
%��gamma�������ж���ÿһ�����ݵ����
    if (cl(ordgamma(i))>0||any(ordgamma(i)==knn(nneigh(ordgamma(i)),:)))
%��ȥ�������ĵ��Լ����ܶȷ�ֵ��(����i����ָ����k�����ڣ�
        continue;
    else
        if (connect(ordgamma(i),nneigh(ordgamma(i)))==1&&cl(nneigh(ordgamma(i)))>0)
%����ָ�����ͨ������ָ���Ҳ�ѷ�����ܶȷ�ֵ�����Aq
            cl(ordgamma(i))=cl(nneigh(ordgamma(i)));
            Aq=[Aq ordgamma(i)];
        else
%����ָ��㲻��ͨ��������ָ�����ͨ��ָ��㲢δ������ܶȷ�ֵ�����Bq
            Bq=[Bq ordgamma(i)];
        end
    end
end

%��Bq��������Aq������ͨ�ĵ����μ���Aq
while ~isempty(Bq)  %��Bq����Ϊ��ʱ
    [i,j,Bq_j,max_d]=mindist(Aq,Bq,dist);
%���Ӻ����ҵ�Aq�����ڵ�i��Bq�����ڵ�j��ʹ����������С��
%�����Bq�е������λ��Bq_j,�Լ�Aq��Bq�����������
    if connect(i,j)==1  %��Aq��Bq�о����������������ͨ��
        cl(j)=cl(i);
        Aq=[Aq j];
        Bq(:,[Bq_j])=[];%ɾ��Bq�еĵ�
    else  %��i���j����ͨ����������Ϊ������
        dist(i,j)=max_d;
    end
    if isempty(Bq)
        break;
    end
    [logic]=link(Aq,Bq,connect);
    if logic==0
        %��Aq�����еĵ���Bq�����ڸ������������ͨ��ʱ���˳�ѭ��
        break;
    end
end
        
%�������������Bq��Ϊ�գ���ʱBq�ڸ��������Aq�ڵĵ������ͨ��
while ~isempty(Bq)
%�ҵ�Aq�����ڵ�i��Bq�����ڵ�j��ʹ����������С�������Bq�е������λ��Bq_j
    [i,j,Bq_j,~]=mindist(Aq,Bq,dist);
    cl(j)=cl(i);
    Bq(:,[Bq_j])=[];%ɾ��Bq�еĵ�
    Aq=[Aq j];%���õ����Aqĩβ
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
%����Ϊֹ����ɾ������   

%���Ϊ��ά���ݼ������ƾ�����
if d==2
cmap=colormap;%���ص�ǰͼ������ɫͼ��ΪRGB��Ԫ��ɵ�3�о���,Ϊ����ɫ��
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

%��һ����ÿ�� cluster
for i=1:NCLUST
  nc=0;  %% �����ۼƵ�ǰ cluster �����ݵ�ĸ���
  for j=1:N
    if (cl(j)==i)
      nc=nc+1;
    end
  end
  fprintf('CLUSTER: %i CENTER: %i ELEMENTS: %i \n', i,icl(i),nc);
end

%�洢������
faa=fopen('CLUSTER_ASSIGNATION','w');
disp('Generated file:CLUSTER_ASSIGNATION')
disp('column 1:element id')
disp('column 2:cluster assignation without noise control')
for i=1:N
   fprintf(faa, '%i %i\n',i,cl(i));
end

%% �������ʶ��
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

%ѡ��һ������Χס�����ľ���
disp('Select a rectangle enclosing cluster noises')

%������ν�ġ�����ͼ2��
%��Ϊָ��һ��λ��
figure('Position',[6 72 scrsz(3)/2 scrsz(4)/1.3]);
subplot(1,1,1)
tt=plot(1:N,theta_sorted(:),'o','MarkerSize',3,'MarkerFaceColor','k','MarkerEdgeColor','k');
%o��ԲȦ����ǵ��С��3����ǵ��ڲ������ɫ����ɫ����ǵ�߿���ɫ����ɫ��
title ('Decision Graph 2','FontSize',10.0)%��ͼ�е�text����
xlabel ('IND')
ylabel ('Sorted \theta')

fig=subplot(1,1,1);
rect = getrect(fig);
%getrect ��ͼ��������ȡһ���������� rect �д�ŵ���
%�������½ǵ����� (x,y) �Լ����ؾ��εĿ�Ⱥ͸߶�

thetamin=rect(2);

%ͳ�������������������
N_NOISE=0;
for i=1:N
  if ( theta(i)>thetamin)
     N_NOISE=N_NOISE+1;
  end
end

%��noise�����ڲ��洢���ı��
for i=1:N_NOISE
    noise(i)=ordtheta(i);
end

%halo(i)=0����õ�Ϊ������
for i=1:N
    halo(i)=cl(i);
    for j=1:N_NOISE
    if i==noise(j)
        halo(i)=0;
    end
    end
end
%����noise�����е�˳���Ա��ں���Ƚ�
noise=find(halo==0);

%����k���ڣ�ʹ��k�����Ż�����ʶ��
%��ʼ��nrho(i)=0,��������ܶ�
nrho=zeros(1,N);
%nrho(i)>0��ʾ��i���ں��ĵ㣬nrho(i)<=0��ʾ��i�������
for i=1:N
    for j=1:k
            if halo(knn(i,j))>0
                nrho(i)=nrho(i)+1;
            else
                nrho(i)=nrho(i)-1;
            end
    end
end
%�洢�п��ܳ�Ϊ���ĵ㣨ȡnrho<=k-1�����ݵ㣩
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
%�ж��Ƿ�Ϊ���,nrho(i)<=0��ʾ��i�������
for i=1:N
    if nrho(i)<=0
        halo(i)=0;%��iΪ���
    else
        halo(i)=cl(i);%��i��Ϊ���
    end
end

NOISE=find(halo==0);

%�ж�NOISE�����Ƿ����仯����δ�����仯��˵����㲻�ڱ仯���˳�ѭ��
if N_NOISE==size(NOISE,2)
    if noise==NOISE
        break;
    end
end
noise=NOISE;
end
fprintf('NUMBER OF NOISES: %i \n', N_NOISE);

if d==2
%����������������ľ�����
cmap=colormap;%���ص�ǰͼ������ɫͼ��ΪRGB��Ԫ��ɵ�3�о���,Ϊ����ɫ��
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

%��һ����ÿ�� cluster
for i=1:NCLUST
  nc=0;  %% �����ۼƵ�ǰ cluster �����ݵ�ĸ���
  nh=0;  %% �����ۼƵ�ǰ cluster �к������ݵ�ĸ���
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

%�洢������
faa=fopen('CLUSTER_ASSIGNATION_NOISE','w');
disp('Generated file:CLUSTER_ASSIGNATION_NOISE')
disp('column 1:element id')
disp('column 2:cluster assignation without noise control')
disp('column 3:cluster assignation with noise control')
for i=1:N
   fprintf(faa, '%i %i %i\n',i,cl(i),halo(i));
end
end
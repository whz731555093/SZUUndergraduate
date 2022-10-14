%% 定位方法1：取三条，分别拟合再平均
clear all;clc;
bx1=0; by1=0;  
bx1=-5e-9; by1=-5e-9;    
bx2=5e-9; by2=5e-9;     
cx=3e-10; cy=3e-10;
d=4e-9;
N=250;
dx0=0.096e-9; dy0=0.096e-9;

xx1 = (-N/2:N/2-1) * dx0; % Input f0 is in natural order
yy1 = (-N/2:N/2-1) * dy0;
[x0,y0]=meshgrid(xx1,yy1);
%y1=exp(-((x0-bx1).^2/(2*cx^2)+(y0-by1).^2/(2*cy^2))) + exp(-((x0-bx2).^2/(2*cx^2)+(y0-by2).^2/(2*cy^2)));
y1=1/(2*pi*cx*cy)*exp(-((x0-bx1).^2/(2*cx^2)+(y0-by1).^2/(2*cy^2))) + 1/(2*pi*cx*cy)*exp(-((x0-bx2).^2/(2*cx^2)+(y0-by2).^2/(2*cy^2)));
figure(1);imshow(y1);

y2=zeros(N,N);
y2(:,N/2)=exp(-xx1/d);
f=conv2(y1,y2);
tj=1;dd=zeros(4,tj);aaa=1;bbb=1;ccc=aaa+bbb+1;dnn=1;
%imshow(mat2gray(f))

A1 = imnoise(mat2gray(f),'poisson');
A2 = imnoise(mat2gray(A1),'gaussian');
A3=wiener2(A2,[3 3]); %step1,filter
%A3 = SF3(1:N,N/2:N+N/2-1);
imshow(A3);
%% cftool拟合x向高斯函数
imagex = sum(A3,1);   %将图像矩阵向着水平方向（x轴）投影叠加
res = double(imagex);   % 将数据转换为double类型
[cols, rows] = size(res);   % 丈量数据矩阵尺寸，此处数据集为1行120列
gsx_axis = 1:size(res,2);   % X轴的坐标区间
gsy_axis = res/max(res);     % Y轴为数据点的值
cftool
%% cftool拟合y向高斯函数
imagey = sum(y1,2);   %将图像矩阵向着水平方向（x轴）投影叠加
res = double(imagey);   % 将数据转换为double类型
[cols, rows] = size(res);   % 丈量数据矩阵尺寸，此处数据集为1行120列
gsx_axis = 1:size(res,1);   % X轴的坐标区间
gsy_axis = res/max(res);     % Y轴为数据点的值
cftool
%% lsqcurvefit拟合高斯函数卷积寿命衰减函数求时间常数的初值
y1=73.92; y2=178.1;
sigmay1=4.419; sigmay2=4.419;
x1=198.0; x2=302.0;
sigmax1=4.118; sigmax2=4.237;
a=round(x1);
b=round(x2);
imagee=A3(:,a)./6+A3(:,a-1)./6+A3(:,a+1)./6+A3(:,b)./6+A3(:,b-1)./6+A3(:,b+1)./6;   %取第最大列进行拟合定时间常数初值
res = double(imagee);   % 将数据转换为double类型
res_y = res';   %将列向量转换为行向量
[rows, cols] = size(res_y);   % 丈量数据矩阵尺寸，此处数据集为1行239列  rows=1，cols=239
ex_axis = 1:size(res_y,2);   %x轴坐标区间
ey_axis = res_y/max(res_y);   %y轴数据归一化

oy(1)=y1;   %高斯函数均值为y1=[52.2527103069531]
oy(2)=y2;   %高斯函数均值为y2=[61.0450913412387]
alocation = find(res_y==max(res_y));   %储存最大值（多个最大值同时存在的情况）及其位置
halfpeaks=find(abs(res_y(:,1:alocation(1))-max(res_y)./2)==min(abs(res_y(:,1:alocation(1))-max(res_y)./2)));  %找到左边的纯高斯的半峰值位置=57
elocation=find(abs(res_y(:,alocation(1):cols)-max(res_y)./2.718281)==min(abs(res_y(:,alocation(1):cols)-max(res_y)./2.718281)));  %%储存时间常数（多个同时存在的情况）及其位置
oy(3)=elocation(1);  %荧光寿命常数初值=[13.4586529387669]
oy(4)=abs(alocation(1)-halfpeaks(1))*2/2.35482;   %标准差=[8.49321816529501];cftool拟合的标准差sigmay1=2.249
oy(5)=0;   %噪声=[0.156884137631035]
oy(6)=1;   %y1振幅初值=[5.85843486636557]
oy(7)=1;   %y2振幅初值=[-3.41001445017942]
fun2=@(oy,x)oy(5)+oy(6)*conv(exp(-((x-oy(1)).^2)./(2*oy(3).^2)),exp(-x./oy(3)),'same')+...
    oy(7)*conv(exp(-((x-oy(2)).^2)./(2*oy(3).^2)),exp(-x./oy(3)),'same');%对高斯函数*单指数衰减函数卷积的函数做拟合
options=optimset('TolX',2.000000e+02); 
[yy,resnormy] = lsqcurvefit(fun2,oy,ex_axis,ey_axis,[],[],options);
yy
%% 曲面拟合
facex_axis = 1:size(A3,2);   %取x的坐标区间（横向120）
facey_axis = 1:size(A3,1);   %取y的坐标区间（纵向239）
[X,Y] = meshgrid(facex_axis,facey_axis);
facez_axis = double(A3);    %将image3数据转换为double类型
xy=[X(:) Y(:)];
z=facez_axis(:);

oi(1)=x1;   %x1=56.8499959303641
oi(2)=x2;   %x2=56.8500000095986
oi(3)=yy(1);   %y1=51.9810100467012
oi(4)=yy(2);   %y2=60.7482279584067
oi(5)=yy(3);   %荧光寿命=13.4586529361055
oi(6)=(sigmax1+sigmay1)/2;   %x1y1的sigma1=2.64248398340200
oi(7)=(sigmax2+sigmay2)/2;  %x2y2的sigma2=3.04349995829075
oi(8)=1;   %振幅=1.02744649527525
oi(9)=1;   %振幅=0.954276775324733
oi(10)=0;   %噪声=0
fun3=@(oi,x) oi(10)+conv(oi(8).*exp(-((x(:,1)-oi(1)).^2+(x(:,2)-oi(3).^2))./(2.*oi(6).^2)),exp(-x(:,2)./oi(5)),'same')+...
    conv(oi(9).*exp(-((x(:,1)-oi(2)).^2+(x(:,2)-oi(4).^2))./(2.*oi(7).^2)),exp(-x(:,2)./oi(5)),'same');
options=optimset('TolX',2.014195e-01); 
[ii,resonrm]=lsqcurvefit(fun3,oi,xy,z,[],[],options);   %resonrm=[5.47346865789654e+179]
ii
%% fitting
sf1=sum(A3);  %A3即image3，sf1即imagex
figure(2);subplot(1,2,1);imshow(A2);subplot(1,2,2);imshow(SF3);
rm=find(sf1==max(sf1));   %找image3的峰值

cm0(1)=xx1(rm);   %？？？将峰值的横坐标储存在cm0（1）
rrm=find(abs(sf1(1:N/2)-max(sf1)/2)==min(abs(sf1(1:N/2)-max(sf1)/2)));   %找到半峰值的位置
cm0(2)=abs(cm0(1)-xx1(rrm))*2/2.35482;%%%标准差
cm0(3)=0;

fun1 = @(cm,x) cm(3)+exp(-(x-cm(1)).^2/(2*cm(2)^2));
options=optimset('TolX',1e-10^5); 
ccm=lsqcurvefit(fun1,[cm0(1) cm0(2) cm0(3)],xx1,sf1,[],[],options);

rmm=find(abs(xx1(1:N/2)-ccm(1))==min(abs(xx1(1:N/2)-ccm(1))));%%rrm啥意思

%寿命衰减拟合
t=1;   cc=zeros(5,ccc);   ee=0;
for i=N/2+rmm-aaa:N/2+rmm+bbb  %%%%%这一步啥意思
Aa=SF3(:,i);   %第i列所有行构成的列向量
A=Aa';  %将第i列向量转为行向量
% figure(3);subplot(1,ccc,t);plot(Aa,'k');hold on;
B=A(1:N);   %由行向量A的第1到第N个行元素组成的子行向量

c0(1)=1;   %start-guess here 
r1=find(B==max(B));   %找到y方向峰值位置记录在r1,B=res_y
c0(2)=xx1(r1(1));  %将峰值的值储存在c0(2),y(1)为峰值的值
rr=find(abs(B(1:r1(1))-max(B)/2)==min(abs(B(1:r1(1))-max(B)/2)));  %找左边半高（类似于找到左边的纯高斯的半峰值）
rrr=find(abs(B(r1(1):N)-max(B)/2.71828)==min(abs(B(r1(1):N)-max(B)/2.71828)));  %右边半高
c0(3)=abs(c0(2)-xx1(rr))*2/2.35482;  %标准差
c0(4)=abs(xx1(rrr+r1(1)-1)-c0(2));  %寿命
c0(5)=0;   %噪声

fun = @(c,x) c(5)+conv(c(1)*exp(-(x-c(2)).^2/(2*c(3)^2)),exp(-x/c(4)));   %对高斯函数*单指数衰减函数卷积的函数做拟合
options=optimset('TolX',1e-10^5); 
[cc,res]=lsqcurvefit(fun,c0,xx1,A,[],[],options);

if (res<1)
    ee=cc+ee;
    t=t+1;%%%这个循环啥意思
end
    Ifit=fun(cc,xx1); %your fitted gaussian in vector 
%     plot(Ifit,'r')
end
% ddcc(dn,1)=res;ddcc(dn,2:4)=cc(2:4);
ddj=ee/(t-1);
if (t>1)
    dd(1:3,dnn)=ddj(2:4);
    dd(4,dnn)=ccm(1);
    dnn=dnn+1;%%%%这个循环啥意思
end

%end


ddd=dd';
xlswrite('dd.xls',ddd);
ds=find(dd(1,:)==0);
if numel(ds)==0%%%查找元素个数
    ds(1)=tj+1;
end
dds(1:(ds(1)-1),1:4)=ddd(1:(ds(1)-1),1:4);
figure()
[locx1,locx2]=hist(dds(:,4),20);%localization accuracy in x
dealfitx(locx2,locx1);
figure()
[locy1,locy2]=hist(dds(:,1),20);%localization accuracy in y
dealfity(locy2,locy1);
% [NN1,XX1]=hist(dd(1,:),10);
figure()
[loct1,loct2]=hist(dds(:,3),20);title('tao= 1000e-9')%localization accuracy in tao
dealfitt(loct2,loct1);
% [NN3,XX3]=hist(dd(3,:),10)
cftool
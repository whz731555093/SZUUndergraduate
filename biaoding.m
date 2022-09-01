clear all;clc;
image0=imread('C:\Users\吴涵哲\Desktop\毕业设计\FLIM标定-2\cal_98\img_000000000__000.tif');
image2=medfilt2(image0);   %进行中值滤波
image3=wiener2(image0,[3 3]);
%% x向坐标标定（水平）
imagex = sum(image3,1);   %将图像矩阵向着水平方向（x轴）投影叠加
res = double(imagex);   % 将数据转换为double类型
[cols, rows] = size(res);   % 丈量数据矩阵尺寸，此处数据集为1行120列
gsx_axis = 1:size(res,2);   % X轴的坐标区间
gsy_axis = res/max(res);     % Y轴为数据点的值
cftool
%% y向坐标标定（竖直）
imagey = sum(image0,2);   %将图像矩阵向着竖直方向（y轴）投影叠加
res = double(imagey);   % 将数据转换为double类型
res_y3 = res';   %将列向量转换为行向量
[cols, rows] = size(res);   % 丈量数据矩阵尺寸，此处数据集为1列239行  rows=1，cols=1040
ex_axis = 1:size(res,1);
ey_axis = res;
%ey_axis = res/max(res);
cftool
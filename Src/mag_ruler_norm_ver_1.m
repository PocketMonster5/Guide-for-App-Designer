%原始数据采集：将位移与8信号数据存入DISP与SIGNAL
load('ExportData_20210429_112112.mat');
DISP=ExportCtrData;SIGNAL=ExportAIData;

%==============================================================================
%标准信号生成
%norm_disp：标准位移，从0按步长step(mm)增至disp_end(mm)
%norm_signal：标准信号，即标准位移相对应的8探测信号，由signal插值且归一化得到
%norm_para：8信号归一化参数
%------------------------------------------------------------------------------
%全局变量：位移步长step(mm)，标准信号数据norm_signal，信号归一化参数norm_para
global step;global norm_signal;global norm_para; %#ok<*GVMIS> 
%标准位移步长与末端点
step=0.01;disp_end=40;
%标准位移数据（未用）
norm_disp=0:step:disp_end;
%遍历SIGNAL生成插值标准信号norm_signal
norm_signal=[];
num=1;%临时变量
for x=0:-step:-disp_end
    while DISP(num)>x
        num=num+1;
    end
    if DISP(num)==x
        norm_signal=[norm_signal;SIGNAL(num,:)];
    elseif DISP(num)<x
        norm_signal=[norm_signal;SIGNAL(num,:)+(SIGNAL(num-1,:)-SIGNAL(num,:))*(x-DISP(num))/(DISP(num-1)-DISP(num))];
    end
end
%幅值归一化
[norm_signal,norm_para]=mapminmax(norm_signal');norm_signal=norm_signal';
%==============================================================================

%==============================================================================
%粗算周期与1，5信号相位差,用于计算迭代初值
%period_n：信号序号的周期；period_disp：位移周期
%phase_1_5：1，5信号相位差，范围0~2*pi，以信号5超前信号1的方向为正
%------------------------------------------------------------------------------
%全局变量：位移周期period_disp(mm)，1、5信号相位差phase_1_5
global period_disp;global phase_1_5;
[~,n_min_1]=min(norm_signal(:,1));[~,n_max_1]=max(norm_signal(:,1));
[~,n_min_5]=min(norm_signal(:,5));[~,n_max_5]=max(norm_signal(:,5));
period_n=2*abs(n_min_1-n_max_1);
period_disp=period_n*step;
if n_min_5>n_min_1
    phase_1_5=(n_min_5-n_min_1)/period_n*2*pi;
else
    phase_1_5=(n_max_5-n_max_1)/period_n*2*pi;
end
%==============================================================================


%==============================================================================
%测试脚本：1000次位移解算
error=zeros(1000,1);
n_test=floor(rand(1000,1)*20000+10000);
for k=1:1000   
    %位移与信号测试数据
    X_true=-ExportCtrData(n_test(k));
    Y_measure=ExportAIData(n_test(k),:);
    Y_measure=Y_measure';
    %调用解算函数
    Disp=Disp_cal(Y_measure);
    %解算误差
    error(k)=Disp-X_true;
end
plot(error);xlabel('测试点');ylabel('位移误差/mm');
m=mean(error)
s=std(error)
%==============================================================================



%% 
%单次位移计算函数
%输入：Y_measure:8*1列向量，即8信号
%输出：Disp:1*1数字，表示距基准数据零点的位移
function [Disp]=Disp_cal(Y_measure)
    %全局变量：位移步长step(mm)，标准信号数据norm_signal，信号归一化参数norm_para
    global step;global norm_signal;global norm_para;
    %全局变量：位移周期period_disp(mm)，1、5信号相位差phase_1_5
    global period_disp;global phase_1_5;

    %信号归一化
    Y_measure=mapminmax('apply',Y_measure,norm_para);

    %计算初值
    %0位移处1，5信号
    y1=norm_signal(1,1);y5=norm_signal(1,5);
    %0位移处1信号的相位
    phi_0=atan2((y5-y1*cos(phase_1_5))/sin(phase_1_5),y1);
    %待测位移处1，5信号
    y1=Y_measure(1);y5=Y_measure(5);
    %待测位移处1信号的相位
    phi_measure=atan2((y5-y1*cos(phase_1_5))/sin(phase_1_5),y1);
    %待测位移的估计值：X
    if phi_measure<phi_0
        phi_measure=phi_measure+2*pi;
    end
    X=period_disp/2/pi*(phi_measure-phi_0);
    
    %高斯牛顿迭代
    step=0.01;%位移步长
    for i=1:10
        num_X=floor(X/step);
        J=(norm_signal(num_X+2,:)-norm_signal(num_X+1,:))/step;
        G=norm_signal(num_X+1,:)+(norm_signal(num_X+2,:)-norm_signal(num_X+1,:))*(X/step-num_X);
        G=G'-Y_measure;J=J';
        X=X+(J'*J)\(J'*(-G));
        %如果迭代过程中位移<0，应补上一个周期长度；
        if X<0
            X=X+period_disp;
            i=1;%由于周期粗略，应重置迭代步数
        end
    end
    Disp=X;
end
# Batch-Effect

# Batch-Effect

### 1.step1-normlize

数据：GSE83521_and_GSE89143
箱线图显示GSE89143第三个样本异常，使用limma::normalizeBetweenArrays处理过，发现效果不如剔除该样本。

![](https://cdn.nlark.com/yuque/0/2019/png/473837/1576486077497-b7603b86-ea1d-464a-b67d-467a06455e85.png#align=left&display=inline&height=1019&name=image.png&originHeight=1019&originWidth=1887&size=122266&status=done&style=none&width=1887)

将该样本剔除。

### 2.使用sva::ComBat()和limma::removeBatchEffect()分别去除批次效应并做对比

前图时sva，后图是limma

#### 1.校正箱线图

![image.png](https://cdn.nlark.com/yuque/0/2019/png/473837/1576487286562-42273ad4-74bc-43c2-932b-4aa86df7e9ce.png#align=left&display=inline&height=570&name=image.png&originHeight=570&originWidth=853&size=60905&status=done&style=none&width=853)

#### 2.PCA图

![image.png](https://cdn.nlark.com/yuque/0/2019/png/473837/1576487335732-7edbc8b1-6b29-42e8-bc8c-5570efb11c6f.png#align=left&display=inline&height=572&name=image.png&originHeight=572&originWidth=849&size=34680&status=done&style=none&width=849)
![image.png](https://cdn.nlark.com/yuque/0/2019/png/473837/1576487478222-f0627a3f-770f-4fd3-9819-bbc50f1014f6.png#align=left&display=inline&height=574&name=image.png&originHeight=574&originWidth=863&size=36006&status=done&style=none&width=863)

#### 3.火山图

![image.png](https://cdn.nlark.com/yuque/0/2019/png/473837/1576487387677-924d43f6-f15c-48be-81d6-146459e9d7a8.png#align=left&display=inline&height=577&name=image.png&originHeight=577&originWidth=856&size=61028&status=done&style=none&width=856)
![image.png](https://cdn.nlark.com/yuque/0/2019/png/473837/1576487516158-b9bc2579-b03e-4978-ac18-5dfd18cd0d27.png#align=left&display=inline&height=579&name=image.png&originHeight=579&originWidth=855&size=61571&status=done&style=none&width=855)

#### 4.差异基因热图

![image.png](https://cdn.nlark.com/yuque/0/2019/png/473837/1576487441582-e31cd262-a86f-465d-b7ef-420bd04ade7e.png#align=left&display=inline&height=1009&name=image.png&originHeight=1009&originWidth=1912&size=181935&status=done&style=none&width=1912)
![image.png](https://cdn.nlark.com/yuque/0/2019/png/473837/1576487554611-89dd67ca-3ab3-46dc-9c04-7403d901ae91.png#align=left&display=inline&height=1015&name=image.png&originHeight=1015&originWidth=1906&size=187539&status=done&style=none&width=1906)

差异基因数量：
sva： 
  down stable     up  
   30   2108     52 
limma：
 down stable     up  
   30   2099     61 
阈值：p<0.05，|logfc|>1
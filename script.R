################################################################################
# パッケージ
################################################################################

# pacman パッケージがインストールされていることを確認する
if (!require("pacman")) install.packages("pacman")

pacman::p_load(flextable,lattice,skimr,broom,broom.mixed,lme4,tidyverse,tidylog)

################################################################################
# データの読み込み
################################################################################

# .RDataは以下から取得したもの
# このなかにsleepstudyが入っている
# https://github.com/kishiyamat/tutorial-lme-vwp/tree/gh-pages
load('.RData')
head(sleepstudy)

# Reaction: 反応時間(連続値応答変数)
# Correct: 回答の正解/不正解(2値応答変数)
# Days: 寝不足何日目か(説明変数)
# Subject: 被験者(ID)

# 例:被験者番号308は、寝不足0日目の結果は、反応時間249msで正解

################################################################################
# 分布確認
################################################################################
skimr::skim(sleepstudy)

# Q.全被験者が0~9日目まで毎日実施したの？ 
# A. YES
sleepstudy %>% 
  group_by(Subject) %>% 
  summarise(
    n = n(),
    min_Days = min(Days),
    max_Days = max(Days),
    n_distinct_Days = n_distinct(Days)
  )

################################################################################
# 線形モデル(単回帰)
################################################################################

# 反応時間と寝不足日数
# formula
f <- Reaction ~ Days
# plot
lattice::xyplot(f, sleepstudy,type=c('p','r')) # plot, regression
# model
m <- lm(f, sleepstudy)
summary(m)

################################################################################
# 一般化線形モデル (ロジスティック回帰)
################################################################################

# 正解率と寝不足日数
# formula
f2 <- Correct ~ Days
# plot
lattice::xyplot(f2, sleepstudy,type=c('p','r'))
# model
m2 <- glm(f2, sleepstudy, family=binomial)
summary(m2)

# odds計算
m2 %>% 
  broom::tidy(exponentiate = T, conf.int=T) %>% 
  mutate(across(where(is.numeric), ~round(.x,digits=2))) 

str(m2)

################################################################################
# 「|」 の意味は「ごとに」
################################################################################

# 「| Subject」は　「Subjectごとに」という意味

f3 = Reaction ~ Days | Subject
lattice::xyplot(f3, sleepstudy, type=c("p", "r"),
       subset=Subject %in% c(370, 371,372))

# 被験者ごとにReactionとDaysの関係
m3 = lme4::lmList(f3, sleepstudy)
summary(m3)

################################################################################
# 線形混合モデル
################################################################################

# 全体的な効果を「固定効果」とよぶ
# 個人差などの構造的なノイズを「ランダム効果」とよぶ

# 知りたいのは,Daysの全体的な効果
# ランダム効果(Days効果の個人差)を除外したい 
# ランダム効果(Days効果の個人差)は (Days | Subject)と表す

# ランダム効果を除外した上で、固定効果を見るために以下のようなモデルを作成する
# outcome ~ 固定効果 + ランダム効果

f4 = Reaction ~ Days + (Days | Subject)
m4 = lme4::lmer(f4, sleepstudy)
summary(m4)

# Linear mixed model fit by REML ['lmerMod']
# Formula: Reaction ~ Days + (Days | Subject)
# Data: sleepstudy
# 
# REML criterion at convergence: 1743.6
# 
# Scaled residuals: 
#   Min      1Q  Median      3Q     Max 
# -3.9536 -0.4634  0.0231  0.4634  5.1793 
# 
# Random effects:
#   Groups   Name        Variance Std.Dev. Corr
# Subject  (Intercept) 612.10   24.741       
# Days         35.07    5.922   0.07
# Residual             654.94   25.592       
# Number of obs: 180, groups:  Subject, 18
# 
# Fixed effects:
#   Estimate Std. Error t value
# (Intercept)  251.405      6.825  36.838
# Days          10.467      1.546   6.771
# 
# Correlation of Fixed Effects:
#   (Intr)
# Days -0.138

################################################################################

# 解釈を引用

# さきほどの lmList では個人ごとにモデリングしたためバラバラでしたが、
# 今回のモデルの summary では Fixed effects: と Random effects: に分かれており、
# 全体的な効果は Fixed effects 側が持っています。
# Fixed effects の Days が 10.467 であることから、 1日増えると反応時間が10.467増えることがわかります。
# ただ Random effects の方は本当に 個人ごとの Daysの効果のバラつきを考慮できているのでしょうか。

# lmList で確認した個人ごとのバラつき、 つまり構造的なノイズは Random effects がモデリングしてます。
# 込み入った話になるのですが、 この「バラバラ感」は Random effects の Std.Dev で確認できます。
# (Intercept) が平均 24.740 ばらついていて、 Days の効果が平均 5.922 ばらついていることからわかります。
# 
# つまり結果は以下のように解釈できます。
# まず、Fixed effects の Days が 10.467 であることから、 1日増えると反応時間が10.467増えることがわかります。
# ただし、その増え具合は被験者で平均 5.922 程度ばらつく、 というものです。
# 切片も同様で、基本は 251.405 だけれども 被験者で平均 24.740 程度ばらつくよ、 という解釈です。

################################################################################

# 参考までに、ランダム効果の実態
lme4::ranef(m4)$Subject

# (Intercept)        Days
# 308   2.2585509   9.1989758

# この結果を用いて、308番さんのReactionを推定する際の係数を以下のように求める
# Daysの係数 = 10.467 + 9.1989758
# 切片 = 251.405 + 2.2585509


################################################################################
# 一般化線形混合モデル
################################################################################

f5 = Correct ~ Days + (Days | Subject)
m5 = lme4::glmer(f5, sleepstudy, family=binomial)
summary(m5)

# odds計算
# 参考
# https://stackoverflow.com/questions/26417005/odds-ratio-and-confidence-intervals-from-glmer-output

# 方法1
cc <- confint(m5,parm='beta_')
ctab <- cbind(est=fixef(m5), cc)
rtab <- exp(ctab)
print(rtab, digits=3)

# 方法2
m5 %>% 
  broom.mixed::tidy(conf.int=TRUE,exponentiate=TRUE,effects="fixed",conf.method="profile") %>% 
  mutate(across(where(is.numeric), ~round(.x,digits=3))) 

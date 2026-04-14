# Transcript Giải Thích Chi Tiết – Bài 4.5: MANOVA 1 Nhân Tố (Dữ Liệu Rùa)

---

## 1. ĐỀ BÀI (Bài 4.5 – Trang 83, Sách Thống Kê Nhiều Chiều Chương 4)

Bài 4.5 sử dụng bộ dữ liệu đo lường hình thể của **24 con rùa cái** và **24 con rùa đực** trên ba biến:

| Biến | Ký hiệu | Đơn vị |
|------|---------|--------|
| Chiều dài mai rùa | $X_1$ = `Length` | mm |
| Chiều rộng mai rùa | $X_2$ = `Width` | mm |
| Chiều cao mai rùa | $X_3$ = `Height` | mm |

**Nhân tố phân nhóm:** Giới tính (`Sex`): female / male → g = 2 nhóm

### Yêu cầu:
- **(a)** Kiểm định sự bằng nhau của **hai vectơ trung bình tổng thể** với $\alpha = 0.01$
- **(b)** Tìm **tổ hợp tuyến tính** của các trung bình thành phần có ý nghĩa quan trọng trong việc bác bỏ $H_0$
- **(c)** Xây dựng **khoảng tin cậy đồng thời** cho hiệu các trung bình thành phần; so sánh với Bonferroni

---

## 2. KIẾN THỨC SỬ DỤNG VÀ VỊ TRÍ TRONG SÁCH

### 2.1. Mô hình MANOVA 1 nhân tố (Trang 65–70)

Đây là bài toán **so sánh vectơ trung bình của g ≥ 2 tổng thể**. Mô hình tổng quát (công thức (1), trang 65):

$$X_{lj} = \mu + \tau_l + \varepsilon_{lj}, \quad j = 1, 2, \ldots, n_l,\ l = 1, 2, \ldots, g$$

Với bài 4.5: $g = 2$ (female/male), $p = 3$ biến (Length, Width, Height), $n_1 = n_2 = 24$.

Giả thuyết kiểm định (trang 68):

$$H_0: \mu_1 = \mu_2 = \cdots = \mu_g = \mu \quad \leftrightarrow \quad H_0: \tau_1 = \tau_2 = \cdots = \tau_g = \mathbf{0}$$

$$H_A: \text{tồn tại một vectơ } \tau_l \neq \mathbf{0}$$

### 2.2. Bảng MANOVA và ma trận B, W (Trang 69)

Phân tích tổng các tích chéo (SSP – Sum of Squares and Products):

$$\underbrace{\sum_{l=1}^{g}\sum_{j=1}^{n_l}(\mathbf{x}_{lj}-\bar{\mathbf{x}})(\mathbf{x}_{lj}-\bar{\mathbf{x}})'}_{\mathbf{B}+\mathbf{W}} = \underbrace{\sum_{l=1}^{g} n_l(\bar{\mathbf{x}}_l - \bar{\mathbf{x}})(\bar{\mathbf{x}}_l - \bar{\mathbf{x}})'}_{\mathbf{B}} + \underbrace{\sum_{l=1}^{g}\sum_{j=1}^{n_l}(\mathbf{x}_{lj}-\bar{\mathbf{x}}_l)(\mathbf{x}_{lj}-\bar{\mathbf{x}}_l)'}_{\mathbf{W}}$$

Trong đó:
- **B** (Between): Ma trận biến thiên giữa các nhóm, bậc tự do = $g - 1$
- **W** (Within/Residual): Ma trận biến thiên trong nhóm, bậc tự do = $\sum n_l - g$

### 2.3. Thống kê Wilks' Lambda (Trang 70)

$$\Lambda^* = \frac{|\mathbf{W}|}{|\mathbf{B} + \mathbf{W}|}$$

Giá trị $\Lambda^*$ nhỏ → bác bỏ $H_0$. Khi $g = 2$ (hai nhóm), có thể chuyển sang phân phối Fisher:

$$\frac{1 - \sqrt{\Lambda^*}}{\sqrt{\Lambda^*}} \cdot \frac{\sum n_g - g - 1}{g - 1} \sim F_{p(g-1),\ \sum n_g - g - 1}$$

*(Trang 71, Ví dụ 4.2 minh hoạ)*

### 2.4. Điều kiện MANOVA (Trang 68)

Ba điều kiện cần kiểm tra:
1. Các mẫu **độc lập** với nhau
2. Các tổng thể có **chung ma trận hiệp phương sai** $\Sigma$ (kiểm định Box's M)
3. Mỗi tổng thể có **phân phối chuẩn nhiều chiều** (kiểm định Shapiro-Wilk đa biến)

### 2.5. Khoảng tin cậy đồng thời Hotelling T² và Bonferroni

Trong trường hợp hai mẫu độc lập, khoảng tin cậy đồng thời $(1-\alpha)$ cho hiệu $\mu_{1i} - \mu_{2i}$ (thành phần thứ $i$):

**T² (simultaneous):**
$$(\bar{x}_{1i} - \bar{x}_{2i}) \pm \sqrt{c^2} \cdot \sqrt{s_{p,ii}\left(\frac{1}{n_1}+\frac{1}{n_2}\right)}$$

Với $c^2 = \frac{(n_1+n_2-2) \cdot p}{n_1+n_2-p-1} \cdot F_{p,\, n_1+n_2-p-1}(1-\alpha)$

**Bonferroni:**
$$(\bar{x}_{1i} - \bar{x}_{2i}) \pm t_{n_1+n_2-2}\!\left(\frac{\alpha}{2p}\right) \cdot \sqrt{s_{p,ii}\left(\frac{1}{n_1}+\frac{1}{n_2}\right)}$$

Trong đó $S_p = \frac{(n_1-1)S_1 + (n_2-1)S_2}{n_1+n_2-2}$ là ma trận hiệp phương sai gộp.

---

## 3. GIẢI THÍCH TỪNG BƯỚC CODE

---

### Bước 0 – Cài đặt và nạp packages

```r
packages <- c("dplyr", "tidyr", "ggplot2", "ggpubr", "rstatix",
              "mvnormtest", "biotools", "car", "knitr")
```

**Lý do từng package:**
- `dplyr`, `tidyr`: thao tác và biến đổi dữ liệu (pivot_longer, filter, group_by)
- `ggplot2`, `ggpubr`: vẽ boxplot và Q-Q plot
- `rstatix`: kiểm định Shapiro-Wilk theo nhóm, Levene test, Tukey HSD
- `mvnormtest`: kiểm định chuẩn **nhiều chiều** (`mshapiro.test`)
- `biotools`: Box's M test (`boxM`) – kiểm tra đồng nhất ma trận hiệp phương sai
- `car`: hàm `Manova()` với nhiều thống kê kiểm định (Wilks, Pillai, ...)

---

### Bước 1 – Nhập và chuẩn bị dữ liệu

```r
turtle <- read.table("exo4-5.dat",
                     header = FALSE,
                     col.names = c("Length", "Width", "Height", "Sex"))

turtle$Sex <- factor(turtle$Sex,
                     levels = c("female", "male"),
                     labels = c("female", "male"))
```

**Giải thích:**
- `read.table()` đọc file `.dat` không có dòng header → đặt tên cột thủ công.
- Biến `Sex` được chuyển thành `factor` (biến phân loại) vì trong MANOVA, nhân tố phân nhóm phải là factor. Nếu để kiểu `character`, hàm `manova()` vẫn hoạt động nhưng dễ gây nhầm lẫn.
- Dữ liệu: 48 quan trắc ($n_1 = n_2 = 24$), $p = 3$ biến liên tục.

---

### Bước 2 – Thống kê mô tả

```r
mo_ta <- turtle %>%
  pivot_longer(cols = c(Length, Width, Height),
               names_to  = "variable",
               values_to = "value") %>%
  group_by(Sex, variable) %>%
  summarise(n = n(), mean = round(mean(value), 3), sd = round(sd(value), 3))
```

**Giải thích:**
- `pivot_longer()` chuyển dữ liệu từ dạng rộng (3 cột biến) sang dạng dài (1 cột `variable`, 1 cột `value`), thuận lợi cho tính thống kê theo nhóm.
- Tính **trung bình mẫu** $\bar{x}_l$ và **độ lệch chuẩn** cho từng biến theo từng giới tính → là ước lượng cho $\mu_l$ trong mô hình (trang 66).
- Đây là bước EDA (Exploratory Data Analysis) nhằm quan sát sơ bộ sự khác biệt giữa female và male trước khi kiểm định.

---

### Bước 3 – Boxplot

```r
ggplot(turtle_long, aes(x = Sex, y = Value, fill = Sex)) +
  geom_boxplot(alpha = 0.7) +
  facet_wrap(~ Variable, scales = "free_y")
```

**Giải thích:**
- `facet_wrap(~ Variable)` tạo 3 ô con riêng biệt cho Length, Width, Height.
- `scales = "free_y"` cho phép trục y có thang đo khác nhau ở mỗi biến (vì đơn vị tương tự nhưng giá trị khác nhau).
- Mục tiêu: quan sát trực quan xem có sự khác biệt rõ ràng giữa female và male không → phù hợp với nhận định trang 73 (Ví dụ 4.4 – Iris): "Ta thấy các chiều dài trung bình ở mỗi loài khác nhau".

---

### Bước 4 – Kiểm tra tính chuẩn từng biến (Shapiro-Wilk đơn biến)

```r
sw_test <- turtle %>%
  group_by(Sex) %>%
  shapiro_test(Length, Width, Height)
```

**Kiến thức:** Điều kiện 3 của MANOVA (trang 68) yêu cầu phân phối chuẩn nhiều chiều. Trước tiên, ta kiểm tra từng biến riêng lẻ:

$$H_0: \text{Dữ liệu có phân phối chuẩn} \quad \text{vs} \quad H_A: \text{Không có phân phối chuẩn}$$

- Kiểm định **Shapiro-Wilk** kiểm tra tính chuẩn cho mẫu nhỏ ($n < 50$).
- Nếu **p-giá trị > 0.05**: không bác bỏ $H_0$ → chấp nhận tính chuẩn.
- Thực hiện theo từng nhóm giới tính (`group_by(Sex)`) vì MANOVA giả định phân phối chuẩn trong **nội bộ từng nhóm**.

---

### Bước 5 – Q-Q Plot

```r
ggqqplot(turtle, x = "Length", facet.by = "Sex",
         title = "Q-Q plot cho biến Length theo giới tính")
```

**Giải thích:**
- Q-Q plot (Quantile-Quantile plot) vẽ các phân vị quan sát so với phân vị lý thuyết của phân phối chuẩn.
- Nếu các điểm **nằm gần đường thẳng tham chiếu**: dữ liệu tuân theo phân phối chuẩn.
- `facet.by = "Sex"` tạo hai ô riêng cho female và male.
- Lặp lại cho cả 3 biến Length, Width, Height.
- Đây là bước **kiểm tra trực quan** bổ sung cho kiểm định Shapiro-Wilk *(tương tự trang 74–75, Ví dụ 4.4)*.

---

### Bước 6 – Kiểm tra chuẩn nhiều chiều (Multivariate Normality)

```r
female_data <- turtle %>% filter(Sex == "female") %>%
  dplyr::select(Length, Width, Height) %>% t()

mshapiro.test(female_data)
```

**Kiến thức:** Điều kiện 3 (trang 68) yêu cầu mỗi tổng thể có phân phối **chuẩn nhiều chiều** $\mathcal{N}_p(\mu_l, \Sigma)$.

- `mshapiro.test()` từ package `mvnormtest` thực hiện kiểm định **Shapiro-Wilk nhiều chiều** (Mardia test) cho toàn bộ vectơ $(X_1, X_2, X_3)$.
- Hàm yêu cầu input là ma trận **p × n** (biến theo hàng, quan trắc theo cột) → cần `t()` (transpose).
- Thực hiện riêng cho `female_data` và `male_data`.
- Nếu **p-giá trị > 0.05**: chấp nhận phân phối chuẩn nhiều chiều *(trang 75: "Do p-giá trị > 0.05, nên giả định dữ liệu có phân phối chuẩn nhiều chiều thỏa")*.

---

### Bước 7 – Kiểm tra đồng nhất ma trận hiệp phương sai: Box's M Test

```r
boxM_result <- boxM(turtle[, c("Length", "Width", "Height")], turtle$Sex)
```

**Kiến thức:** Điều kiện 2 của MANOVA (trang 68): các tổng thể có **chung ma trận hiệp phương sai** $\Sigma$.

Kiểm định giả thuyết (Phần 4.2, trang 72):

$$H_0: \Sigma_1 = \Sigma_2 \quad \text{vs} \quad H_A: \Sigma_1 \neq \Sigma_2$$

- `boxM()` từ package `biotools` tính thống kê Box's M:

$$M = \left[\sum_l (n_l - 1)\right]\ln|S_{\text{pooled}}| - \sum_l (n_l-1)\ln|S_l|$$

Sau đó hiệu chỉnh thành $C = (1-u) \cdot M$ có phân phối xấp xỉ $\chi^2$ *(công thức trang 73)*.

- Nếu **p-giá trị > 0.05**: chấp nhận đồng nhất ma trận hiệp phương sai → MANOVA hợp lệ.
- **Lưu ý thực tế:** Box's M rất nhạy cảm với cỡ mẫu lớn. Trang 75 ghi rõ: "vì ta thực hiện MANOVA 1 nhân tố trong trường hợp cân bằng ($n_1 = n_2 = n_3 = 50$) nên ta vẫn có thể tiếp tục quá trình phân tích".

---

### Bước 8 – Kiểm tra đồng nhất phương sai từng biến: Levene Test

```r
levene_test(turtle, Length ~ Sex)
levene_test(turtle, Width  ~ Sex)
levene_test(turtle, Height ~ Sex)
```

**Giải thích:**
- Levene test kiểm định $H_0: \sigma^2_{\text{female}} = \sigma^2_{\text{male}}$ cho **từng biến riêng lẻ** (phiên bản đơn chiều của Box's M).
- Đây là bước bổ sung để xem biến nào có phương sai không đồng nhất.
- Nếu Levene test cho thấy vi phạm nhưng Box's M không có ý nghĩa thống kê mạnh, và mẫu cân bằng → MANOVA vẫn bền vững.

---

### Bước 9 (Câu a) – Kiểm định MANOVA bằng Wilks' Lambda

```r
model <- lm(cbind(Length, Width, Height) ~ Sex, data = turtle)
summary(manova(model), test = "Wilks")
```

**Kiến thức cốt lõi (Trang 70):**

**Bước 9.1 – Xây dựng mô hình hồi quy đa biến:**
- `lm(cbind(Length, Width, Height) ~ Sex)`: hồi quy tuyến tính với biến phụ thuộc là **ma trận** gồm 3 cột (p = 3 biến). R tự động tính $\bar{\mathbf{x}}_l$ và ma trận B, W.

**Bước 9.2 – Kiểm định MANOVA:**
- `manova(model)` tính thống kê Wilks' Lambda:

$$\Lambda^* = \frac{|\mathbf{W}|}{|\mathbf{B} + \mathbf{W}|}$$

- `test = "Wilks"` chỉ định dùng thống kê Wilks' Lambda (mặc định, phù hợp với lý thuyết trang 70).
- Với $g = 2$ nhóm, $p = 3$ biến: R tự động chuyển $\Lambda^*$ sang thống kê F xấp xỉ:

$$F = \frac{1 - \sqrt{\Lambda^*}}{\sqrt{\Lambda^*}} \cdot \frac{n - g - 1}{g - 1} \sim F_{p(g-1),\ n-g-1}$$

*(Công thức trang 71, Ví dụ 4.2)*

**Kết luận:** Nếu **p-giá trị < 0.01**: bác bỏ $H_0$, kết luận hai vectơ trung bình tổng thể của rùa cái và rùa đực **khác nhau có ý nghĩa thống kê**.

---

### Bước 10 – Post-hoc: ANOVA từng biến

```r
summary.aov(model)
```

**Giải thích:**
- Sau khi MANOVA bác bỏ $H_0$, ta muốn biết **biến nào** đóng góp vào sự khác biệt.
- `summary.aov(model)` thực hiện 3 kiểm định ANOVA đơn biến riêng lẻ:
  - $F$-test cho `Length ~ Sex`
  - $F$-test cho `Width ~ Sex`
  - $F$-test cho `Height ~ Sex`
- Đây là phiên bản đơn chiều của MANOVA (phần 4.1, trang 65–68): mỗi kiểm định dùng thống kê:

$$F = \frac{SS_{\text{tr}}/(g-1)}{SS_{\text{res}}/(\sum n_l - g)} \sim F_{g-1,\ \sum n_l - g}$$

- **Lưu ý:** ANOVA từng biến có lạm phát lỗi loại I nếu không điều chỉnh; so sánh với kết quả MANOVA để kết luận chắc chắn hơn.

---

### Bước 11 – Post-hoc: So sánh theo cặp Tukey HSD

```r
fit1 <- aov(Length ~ Sex, data = dat_clean)
TukeyHSD(fit1)
```

**Giải thích:**
- Với g = 2 nhóm (female vs male), so sánh cặp Tukey chỉ có 1 cặp duy nhất → kết quả tương đương ANOVA đơn biến.
- Tukey HSD (Honestly Significant Difference) điều chỉnh p-giá trị để kiểm soát **Family-Wise Error Rate** khi so sánh nhiều cặp.
- Lặp lại cho Width và Height.
- Nếu **p-adj < 0.01**: biến đó có sự khác biệt có ý nghĩa thống kê giữa female và male.

---

### Bước 12 (Câu b) – Tổ hợp tuyến tính quan trọng: Phân tích Canonical Discriminant

```r
library(candisc)
fit_b <- lm(cbind(Length, Width, Height) ~ Sex, data = dat_b)
can   <- candisc(fit_b)
print(coef(can, type = "raw"))
```

**Kiến thức (Trang 68–70):**

Bài toán tìm tổ hợp tuyến tính $\mathbf{a'x}$ sao cho **sự phân biệt giữa hai nhóm là lớn nhất**. Đây tương đương với tìm vectơ $\mathbf{a}$ tối đa hóa:

$$\frac{\mathbf{a'Ba}}{\mathbf{a'Wa}}$$

Nghiệm của bài toán là **eigenvector** ứng với eigenvalue lớn nhất của ma trận $\mathbf{W}^{-1}\mathbf{B}$.

- `candisc()` (Canonical Discriminant Analysis) thực hiện phép phân tích này.
- `coef(can, type = "raw")`: trả về **hệ số thô** (raw canonical coefficients) $\mathbf{a} = (a_1, a_2, a_3)'$ cho tổ hợp $a_1 \cdot \text{Length} + a_2 \cdot \text{Width} + a_3 \cdot \text{Height}$.
- Biến nào có **|hệ số| lớn hơn** → đóng góp nhiều hơn vào việc phân biệt hai nhóm → quan trọng hơn trong việc bác bỏ $H_0$.

---

### Bước 13 (Câu c) – Khoảng tin cậy đồng thời T² và Bonferroni

```r
# Tách dữ liệu theo nhóm
x1 <- dat_c[dat_c$Sex == "male",   c("Length", "Width", "Height")]
x2 <- dat_c[dat_c$Sex == "female", c("Length", "Width", "Height")]

n1 <- nrow(x1); n2 <- nrow(x2); p <- ncol(x1)

xbar1 <- colMeans(x1)
xbar2 <- colMeans(x2)

S1 <- cov(x1); S2 <- cov(x2)
Sp <- ((n1 - 1) * S1 + (n2 - 1) * S2) / (n1 + n2 - 2)
```

**Kiến thức:**
- `colMeans()` tính vectơ trung bình mẫu $\bar{\mathbf{x}}_1, \bar{\mathbf{x}}_2$.
- `cov()` tính ma trận hiệp phương sai mẫu $S_1, S_2$.
- `Sp` là **ma trận hiệp phương sai gộp** (pooled covariance matrix):

$$S_p = \frac{(n_1-1)S_1 + (n_2-1)S_2}{n_1+n_2-2}$$

Đây là ước lượng cho ma trận $\Sigma$ chung (điều kiện 2 của MANOVA, trang 68).

---

```r
d  <- xbar1 - xbar2
se <- sqrt(diag(Sp) * (1/n1 + 1/n2))
```

- `d`: vectơ **hiệu trung bình mẫu** $\bar{x}_{1i} - \bar{x}_{2i}$ cho từng biến $i$.
- `diag(Sp)`: lấy các **phần tử đường chéo** của $S_p$ (phương sai gộp của từng biến).
- `se`: **sai số chuẩn** của hiệu trung bình mỗi biến:

$$\text{se}_i = \sqrt{s_{p,ii} \cdot \left(\frac{1}{n_1} + \frac{1}{n_2}\right)}$$

---

**Khoảng tin cậy đồng thời T²:**

```r
c2     <- ((n1 + n2 - 2) * p / (n1 + n2 - p - 1)) * qf(1 - alpha, p, n1 + n2 - p - 1)
lwr_t2 <- d - sqrt(c2) * se
upr_t2 <- d + sqrt(c2) * se
```

**Công thức:**

$$c^2 = \frac{(n_1+n_2-2) \cdot p}{n_1+n_2-p-1} \cdot F_{p,\ n_1+n_2-p-1}(1-\alpha)$$

- `qf(1 - alpha, p, n1 + n2 - p - 1)`: phân vị $1-\alpha$ của phân phối Fisher với bậc tự do $p$ và $n_1+n_2-p-1$.
- Với $\alpha = 0.01$, $p = 3$, $n_1+n_2-p-1 = 48-3-1 = 44$: tính $F_{3,44}(0.99)$.
- Khoảng tin cậy cho thành phần thứ $i$:

$$\left[(\bar{x}_{1i}-\bar{x}_{2i}) - \sqrt{c^2} \cdot \text{se}_i,\ (\bar{x}_{1i}-\bar{x}_{2i}) + \sqrt{c^2} \cdot \text{se}_i\right]$$

Nếu **0 không nằm trong khoảng**: thành phần thứ $i$ có sự khác biệt có ý nghĩa.

---

**Khoảng tin cậy Bonferroni:**

```r
t_bon   <- qt(1 - alpha / (2 * p), df = n1 + n2 - 2)
lwr_bon <- d - t_bon * se
upr_bon <- d + t_bon * se
```

**Công thức:**

$$t_{n_1+n_2-2}\!\left(\frac{\alpha}{2p}\right)$$

- `qt(1 - alpha/(2*p), df = n1+n2-2)`: phân vị $1 - \alpha/(2p)$ của phân phối t với $n_1+n_2-2 = 46$ bậc tự do.
- Bonferroni chia mức ý nghĩa $\alpha$ đều cho $p = 3$ biến → $\alpha/p = 0.01/3$ cho mỗi biến → kiểm soát **Family-Wise Error Rate**.
- Khoảng Bonferroni thường **hẹp hơn** khoảng T² vì dùng phân phối t (1 bậc tự do tốt hơn) thay vì F, đặc biệt khi số biến $p$ nhỏ.

---

**So sánh hai phương pháp:**

| Phương pháp | Ngưỡng | Độ rộng khoảng |
|-------------|--------|----------------|
| T² (đồng thời) | $F_{p,n-p-1}$ | Rộng hơn (bảo thủ) |
| Bonferroni | $t_{n-2}(\alpha/2p)$ | Hẹp hơn (chặt hơn) |

Kết quả kỳ vọng từ code: **cả ba khoảng cho Length, Width, Height đều không chứa 0** → khẳng định rùa cái và rùa đực có kích thước mai rùa khác nhau ở cả 3 chiều đo, với mức ý nghĩa 1%.

---

## 4. TỔNG KẾT QUY TRÌNH PHÂN TÍCH

```
Nhập dữ liệu
    ↓
Thống kê mô tả + Boxplot
    ↓
Kiểm tra điều kiện MANOVA (Trang 68)
    ├── Shapiro-Wilk đơn biến (từng nhóm)
    ├── Q-Q Plot (trực quan)
    ├── Shapiro-Wilk nhiều chiều (mshapiro.test)
    └── Box's M Test + Levene Test (đồng nhất Σ)
    ↓
MANOVA – Kiểm định Wilks' Lambda (Trang 70)
    ↓ [Bác bỏ H₀]
Post-hoc phân tích
    ├── ANOVA từng biến (summary.aov)
    ├── Tukey HSD so sánh cặp
    ├── Canonical Discriminant (tổ hợp tuyến tính tối ưu)
    └── KTC đồng thời T² và Bonferroni
```

---

## 5. THAM CHIẾU SÁCH

| Nội dung | Trang |
|----------|-------|
| Mô hình MANOVA 1 nhân tố, mô hình tổng quát | 65–66 |
| Phân tích tổng bình phương SSP, bảng MANOVA | 67–69 |
| Thống kê Wilks' Lambda, bác bỏ H₀ | 70 |
| Ví dụ 4.2 – công thức F từ Wilks' Lambda | 71 |
| Điều kiện MANOVA (độc lập, cùng Σ, chuẩn nhiều chiều) | 68 |
| Kiểm định Box's M, phần 4.2 | 72–73 |
| Ví dụ 4.4 – MANOVA với dữ liệu Iris (minh hoạ đầy đủ) | 73–76 |
| Khoảng tin cậy đồng thời T² và Bonferroni | Chương 3 (áp dụng cho hai mẫu) |

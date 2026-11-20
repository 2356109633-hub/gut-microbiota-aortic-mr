# GitHub + Zenodo: 给代码仓库拿 DOI

目标：把 `Code_repository` 放到 GitHub，打一个 `v1.0.0` release，让 Zenodo 自动归档并生成 DOI，然后把 DOI 填回论文。

> 不要直接把 GitHub 网页链接当 DOI 用。BMC 的 Data availability 里写 `https://doi.org/10.5281/zenodo.XXXXXXX` 最稳。

---

## 方案 A：GitHub + Zenodo（推荐）

### 1. 准备 GitHub 账号和 Git

- GitHub：https://github.com
- Git for Windows：https://git-scm.com/download/win

装完 Git 后，在 PowerShell 里确认：

```powershell
git --version
```

### 2. 在 GitHub 建一个空仓库

登录 GitHub → 右上角 **+** → **New repository**

建议：

- Repository name：`gut-microbiota-aortic-mr`
- Visibility：**Public**
- 不要勾选 README、.gitignore、license（本地已经有了）

建完后记下仓库地址，例如：

```text
https://github.com/<你的用户名>/gut-microbiota-aortic-mr.git
```

### 3. 本地初始化并推送

在 PowerShell 里：

```powershell
cd "C:\Users\Hall of fame\Desktop\mr-PDR数据\BMC Microbiology\Code_repository"

git init
git add .
git commit -m "Code for BMC Microbiology submission v1.0.0"
git branch -M main
git remote add origin https://github.com/<你的用户名>/gut-microbiota-aortic-mr.git
git push -u origin main
```

如果第一次 push 要求登录，按 GitHub 提示走浏览器授权。

### 4. 确认仓库里有什么

应该包含：

```text
analysis/
data/
results/
.zenodo.json
CITATION.cff
LICENSE
README.md
ZENODO_UPLOAD_GUIDE.md
```

不应该包含：

```text
finngen_R13_*.gz
AAAgen 原始 meta-analysis 文件
任何个体级数据
```

当前 `Code_repository` 已经排除了这些大文件，直接推即可。

### 5. 在 Zenodo 连接 GitHub

1. 打开 https://zenodo.org
2. 右上角 **Sign in** → **Sign in with GitHub**
3. 授权后，右上角头像 → **Settings**
4. 左侧 **GitHub**
5. 找到 `gut-microbiota-aortic-mr`
6. 把开关拨到 **ON**

如果列表里没有仓库，点 **Sync now**。

### 6. 在 GitHub 创建 release

两种方式，任选：

**网页方式：**

1. 打开 GitHub 仓库
2. 右侧 **Releases** → **Draft a new release**
3. Tag：`v1.0.0`
4. Title：`v1.0.0`
5. 点 **Publish release**

**命令行方式：**

```powershell
cd "C:\Users\Hall of fame\Desktop\mr-PDR数据\BMC Microbiology\Code_repository"
git tag -a v1.0.0 -m "v1.0.0"
git push origin v1.0.0
```

### 7. 等 Zenodo 归档

Release 发布后，Zenodo 通常几分钟内自动归档。

回到 Zenodo → **Upload** 或右上角头像 → **My dashboard**，应该能看到：

```text
Gut microbiota and aortic aneurysm/dissection...
Version 1.0.0
DOI: 10.5281/zenodo.XXXXXXX
```

`XXXXXXX` 就是你的 DOI。

### 8. 复制 DOI

格式：

```text
https://doi.org/10.5281/zenodo.XXXXXXX
```

Zenodo 页面通常有两个：

- **DOI**：这一版 `v1.0.0` 的 DOI；
- **Concept DOI**：所有版本的 DOI。

论文里引用 **这一版的 DOI** 即可；以后代码更新，Concept DOI 不变。

### 9. 填回论文和材料

拿到 DOI 后，需要替换这些位置：

```text
[GitHub/Zenodo DOI to be inserted]
```

具体文件：

- 主稿 `Availability of data and materials`
- `Cover_letter_BMC_Microbiology.docx` 第 5 段
- `Code_repository/README.md`
- `Code_repository/CITATION.cff`（可加 `doi:` 字段）
- `BMC_submission_metadata.md`
- `BMC_submission_materials_checklist.md`

替换后重新生成：

- 主稿 DOCX
- Cover letter DOCX + PDF
- 投稿材料包 ZIP

---

## 方案 B：只用 Zenodo（不碰 GitHub）

如果不想装 Git：

1. 打开 https://zenodo.org
2. 注册/登录
3. 点 **New upload**
4. 上传 `Code_repository_v1.0.0.zip`
5. 填：
   - Title：论文标题
   - Authors：10 位作者，按顺序
   - Description：一段仓库说明
   - License：MIT
   - Version：`1.0.0`
   - Keywords：Mendelian randomization; gut microbiota; aortic aneurysm; aortic dissection
6. **Publish**
7. 复制 DOI

缺点：以后更新代码要手动传新版本；GitHub 方式更自动。

---

## 方案 C：GitHub 但不接 Zenodo

只能得到 GitHub 链接，没有 DOI。

BMC 可以接受 GitHub URL，但 Data availability 里有一个正式 DOI 更好，尤其是审稿人点开就能看归档版本。

---

## 常见问题

### 1. Zenodo 没自动归档

检查：

- GitHub 仓库是不是 **Public**；
- Zenodo Settings → GitHub 里仓库开关是不是 ON；
- release 是不是正式 **Published**，不是 draft；
- 等几分钟，再刷新 Zenodo dashboard。

### 2. DOI 有了，但作者/标题不对

在 Zenodo 记录页面点 **Edit**，改 metadata 后保存。已发布的记录可以编辑 metadata，DOI 不变。

### 3. 代码里还有 `[Zenodo DOI to be inserted]`

正常。拿到 DOI 后统一替换，不要只改一处。

### 4. 要不要把 DOI 写进 `CITATION.cff`

可以。在 `CITATION.cff` 里加：

```yaml
doi: "10.5281/zenodo.XXXXXXX"
```

### 5. 以后代码改了怎么办

- 小改：GitHub 新建 release `v1.0.1`，Zenodo 生成新 DOI；
- 论文里引用最初归档的 `v1.0.0` DOI，或引用 Concept DOI；
- 不要覆盖已发布的 Zenodo 版本。

---

## 最短路径

```text
GitHub 建 public repo
→ 本地 git push
→ Zenodo 连 GitHub
→ GitHub 发 v1.0.0 release
→ Zenodo 自动给 DOI
→ 把 DOI 填回论文
```

需要我帮你把 DOI 填回主稿、cover letter 和材料包，拿到 `10.5281/zenodo.XXXXXXX` 后发我即可。

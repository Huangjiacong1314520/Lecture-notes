本地已经安装Git。

宏碁电脑的SSH key已经上传。

# Github上传本地文件流程：

首先在Github上要新建一个库，相当于远程仓库，用来放你这次上传的文件夹以及文件，假设其SSH是git@github.com:Huangjiacong1314520/Lecture-notes.git

比如我想上传A文件夹里的东西，我就在A文件夹里打开命令行，然后

### 如果本地文件夹还不是 Git 仓库：

**初始化 Git 并连接远程仓库**：

```bash
git init
git add .
git commit -m "first commit"
git branch -M main
git remote add origin git@github.com:Huangjiacong1314520/Lecture-notes.git
git push -u origin main
```

### 如果本地文件夹已经是 Git 仓库：

```bash
git remote add origin git@github.com:Huangjiacong1314520/Lecture-notes.git
git branch -M main
git push -u origin main
```

# 当你本地进行修改后，上传到 GitHub 的流程如下：

# 标准上传流程：

### 1. **查看修改状态**
```bash
git status
```
这会显示哪些文件被修改、新增或删除。

### 2. **添加修改到暂存区**
```bash
# 添加所有修改
git add .

# 或者添加特定文件
git add 文件名
```

### 3. **提交修改**
```bash
git commit -m "描述你的修改内容"
```
提交信息应该清晰说明这次修改的目的，例如：
- `"添加运动控制基础章节"`
- `"修复代码示例错误"`
- `"更新图片资源"`

### 4. **推送到远程仓库**
```bash
git push origin main
```

## 完整示例：
```bash
# 1. 查看状态
git status

# 2. 添加所有修改
git add .

# 3. 提交修改
git commit -m "更新运动控制算法文档"

# 4. 推送到GitHub
git push origin main
```

## 其他常用命令：

### 查看提交历史
```bash
git log
```

### 查看具体修改内容
```bash
git diff
```

### 如果只想提交部分文件
```bash
# 逐个选择要提交的文件
git add -p
```

### 撤销未提交的修改
```bash
# 撤销某个文件的修改
git checkout -- 文件名

# 撤销所有未提交的修改
git checkout -- .
```

## 工作流程总结：
1. **修改文件** → 2. **git add** → 3. **git commit** → 4. **git push**

每次完成一些有意义的修改后，就执行这个流程，保持代码的版本管理。

现在你可以尝试修改一些文件，然后按照这个流程上传！
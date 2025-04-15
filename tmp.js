// Fully implemented version: Updated New Chat button logic.
document.addEventListener('DOMContentLoaded', () => {
    // --- Element References & Null Checks ---
    const chatContainer = document.getElementById('chat-container');
    const chatForm = document.getElementById('chat-form');
    const promptInput = document.getElementById('prompt-input');
    const fileInput = document.getElementById('file-input');
    const filePreviewArea = document.getElementById('file-preview-area');
    const sendButton = document.getElementById('send-button');
    const sendIcon = sendButton ? sendButton.querySelector('.fa-arrow-up') : null;
    const loadingIndicator = document.getElementById('loading-indicator');
    const themeToggleButton = document.getElementById('theme-toggle');
    const settingsButton = document.getElementById('settings-button');
    const settingsModalOverlay = document.getElementById('settings-modal-overlay');
    const settingsModalContent = settingsModalOverlay ? settingsModalOverlay.querySelector('.modal-content') : null;
    const settingsModelInput = document.getElementById('setting-model-name');
    const settingsTempSlider = document.getElementById('setting-temperature');
    const settingsTempValueDisplay = document.getElementById('temperature-value-display');
    const settingsMaxTokensInput = document.getElementById('setting-max-tokens');
    const settingsSaveButton = document.getElementById('settings-save-button');
    const settingsCloseButton = document.getElementById('settings-close-button');
    const settingsToolDropdown = document.getElementById('setting-enabled-tool'); // Added Tool Dropdown
    const headerModelDisplay = document.getElementById('model-name-display');
    const footerModelDisplay = document.getElementById('footer-model-name');
    const footerTempDisplay = document.getElementById('footer-temp-display');
    const newChatButton = document.getElementById('new-chat-button'); // New Chat Button
    const passwordModalOverlay = document.getElementById('password-modal-overlay');
    const passwordForm = document.getElementById('password-form');
    const passwordInput = document.getElementById('password-input');
    const passwordSubmitButton = document.getElementById('password-submit-button');
    const passwordError = document.getElementById('password-error');

    // Check if essential elements were found
    const essentialElements = {
        chatContainer, chatForm, promptInput, fileInput, filePreviewArea, sendButton, sendIcon, loadingIndicator,
        themeToggleButton, settingsButton, settingsModalOverlay, settingsModalContent, settingsModelInput,
        settingsTempSlider, settingsTempValueDisplay, settingsMaxTokensInput, settingsSaveButton, settingsCloseButton,
        settingsToolDropdown, // Added Tool Dropdown
        headerModelDisplay, footerModelDisplay, footerTempDisplay, newChatButton, passwordModalOverlay,
        passwordForm, passwordInput, passwordSubmitButton, passwordError
    };

    let allElementsFound = true;
    for (const [name, element] of Object.entries(essentialElements)) {
        if (!element && name !== 'passwordModalOverlay' && name !== 'passwordForm' && name !== 'passwordInput' && name !== 'passwordSubmitButton' && name !== 'passwordError') { // Password elements are optional
            console.error(`Initialization Error: Element with ID/selector for '${name}' not found.`);
            allElementsFound = false;
        }
    }
    if (!chatForm || !promptInput || !chatContainer) {
        if(typeof alert !== 'undefined') alert("Critical UI elements are missing. App cannot start."); return;
    }

    // --- State Variables ---
    let attachedFiles = [];
    let typingIndicatorElement = null;
    const initialModelName = headerModelDisplay?.textContent || "gemini-2.5-pro-preview-03-25";
    let currentSettings = {
        modelName: initialModelName,
        temperature: 0.7,
        maxOutputTokens: 20000,
        enabledTool: 'search' // Default tool
    };
    let isAuthenticated = !window.passwordRequired;

    // --- Configure Marked & Highlight.js ---
    if (typeof marked !== 'undefined') {
        marked.setOptions({
            gfm: true, breaks: true,
            highlight: (typeof hljs !== 'undefined') ? function(code, lang) {
                const language = hljs.getLanguage(lang) ? lang : 'plaintext';
                try { return hljs.highlight(code, { language, ignoreIllegals: true }).value; }
                catch (error) { console.error("Highlight.js error:", error); try { return hljs.highlight(code, { language: 'plaintext', ignoreIllegals: true }).value; } catch (e) { return code.replace(/</g, "&lt;").replace(/>/g, "&gt;");}}
            } : undefined
        });
        if (typeof hljs === 'undefined') { console.warn("highlight.js not loaded."); }
    } else { console.error("Error: marked.js not loaded!"); }

    // --- Initial Setup ---
    loadSettings();
    updateSendButtonState();
    addInitialGreeting();
    scrollToBottom();
    if (window.passwordRequired && !isAuthenticated) {
        if (passwordModalOverlay) passwordModalOverlay.classList.add('active');
        disableChatInput(true);
    } else {
        document.body.classList.remove('password-protected');
        disableChatInput(false);
    }

    // --- Theme Toggle ---
    const applyTheme = (theme) => {
        if (theme === 'dark') { document.documentElement.classList.add('dark'); }
        else { document.documentElement.classList.remove('dark'); }
    };
   const savedTheme = localStorage.getItem('theme') || 'dark'; // Default to dark
   applyTheme(savedTheme);
   if (themeToggleButton) {
       themeToggleButton.addEventListener('click', () => {
            const isDarkMode = document.documentElement.classList.contains('dark');
            const newTheme = isDarkMode ? 'light' : 'dark';
            applyTheme(newTheme);
            localStorage.setItem('theme', newTheme);
        });
   } else { console.error("Theme toggle button not found"); }


   // --- Auto-resize Textarea ---
   if (promptInput) {
       promptInput.addEventListener('input', () => {
           promptInput.style.height = 'auto';
           promptInput.style.height = `${promptInput.scrollHeight}px`;
           if (promptInput.scrollHeight > 200) { // Max height example
               promptInput.style.overflowY = 'auto';
               promptInput.style.height = '200px';
           } else {
               promptInput.style.overflowY = 'hidden';
           }
           updateSendButtonState();
       });
   } else { console.error("Prompt input element not found"); }

    // --- File Handling (Upload Button & Paste) ---
    if (fileInput) { fileInput.addEventListener('change', handleFileSelection); }
    else { console.error("File input element not found"); }
    if (promptInput) { promptInput.addEventListener('paste', handlePaste); }

    function handleFileSelection(event) {
        const files = event.target?.files;
        if (!files) return;
        const MAX_FILES = 5;
        if (attachedFiles.length >= MAX_FILES) { alert(`Max ${MAX_FILES} files.`); return; }
        let filesAdded = 0;
        for (const file of files) {
             if (attachedFiles.length >= MAX_FILES) { alert(`File limit reached.`); break; }
            if (addFileToList(file)) filesAdded++;
        }
        if (filesAdded > 0) updateSendButtonState();
        if (event.target) event.target.value = '';
    }
    function handlePaste(event) {
        const items = event.clipboardData?.items;
        if (!items) return;
        let imageFoundAndAdded = false;
        const MAX_FILES = 5;
        for (const item of items) {
            if (item.kind === 'file' && item.type.startsWith('image/')) {
                 if (attachedFiles.length >= MAX_FILES) { alert(`File limit reached.`); continue; }
                const file = item.getAsFile();
                if (file) { if (addFileToList(file)) imageFoundAndAdded = true; }
                else { console.error("Could not get file from clipboard item."); }
            }
        }
        if (imageFoundAndAdded) { event.preventDefault(); }
    }
     function addFileToList(file) {
        if (!file) return false;
        if (!attachedFiles.some(f => f.name === file.name && f.size === file.size)) {
            attachedFiles.push(file);
            createFilePreview(file);
            updateSendButtonState();
            return true;
        }
        return false;
    }
    function createFilePreview(file) {
        if (!filePreviewArea) { console.error("Preview area not found."); return; }
        const previewElement = document.createElement('div');
        previewElement.className = 'file-preview bg-gray-200 dark:bg-gray-700 p-1 px-2 rounded-lg flex items-center space-x-2 text-xs';
        previewElement.dataset.fileName = file.name;
        let iconClass = 'fa-file';
        if (file.type.startsWith('image/')) iconClass = 'fa-file-image';
        else if (file.type === 'application/pdf') iconClass = 'fa-file-pdf';
        else if (file.type.startsWith('text/')) iconClass = 'fa-file-alt';
        else if (file.type.includes('python') || file.name.endsWith('.py')) iconClass = 'fa-file-code';
        previewElement.innerHTML = `<i class="fas ${iconClass} text-gray-600 dark:text-gray-400"></i> <span class="truncate max-w-[100px]">${file.name}</span> <button type="button" class="remove-file-btn text-red-500 hover:text-red-700 ml-auto text-lg leading-none" title="Remove file">&times;</button>`;
        const removeBtn = previewElement.querySelector('.remove-file-btn');
        if (removeBtn) { removeBtn.addEventListener('click', () => { removeFile(file.name); previewElement.remove(); }); }
        filePreviewArea.appendChild(previewElement);
    }
    function removeFile(fileName) {
        attachedFiles = attachedFiles.filter(f => f.name !== fileName);
        updateSendButtonState();
    }
    function clearFilePreviews() {
        if (filePreviewArea) filePreviewArea.innerHTML = '';
        attachedFiles = [];
        if (fileInput) fileInput.value = '';
        updateSendButtonState();
    }
    function updateSendButtonState() {
        if (!promptInput || !sendButton) return;
        const hasText = promptInput.value.trim().length > 0;
        const hasFiles = attachedFiles.length > 0;
        const isDisabled = !hasText && !hasFiles;
        sendButton.disabled = isDisabled;
        if (isDisabled) { sendButton.classList.remove('bg-gemini-blue', 'text-white', 'hover:bg-blue-600', 'cursor-pointer'); sendButton.classList.add('bg-gray-300', 'dark:bg-gray-500', 'text-gray-700', 'dark:text-gray-200', 'cursor-not-allowed', 'opacity-50'); }
        else { sendButton.classList.remove('bg-gray-300', 'dark:bg-gray-500', 'text-gray-700', 'dark:text-gray-200', 'cursor-not-allowed', 'opacity-50'); sendButton.classList.add('bg-gemini-blue', 'text-white', 'hover:bg-blue-600', 'cursor-pointer'); }
    }

    // --- Settings Modal Logic ---
    function openSettingsModal() {
        if (!settingsModalOverlay || !settingsModelInput || !settingsTempSlider || !settingsTempValueDisplay || !settingsMaxTokensInput || !settingsToolDropdown) {
             console.error("Cannot open settings: Modal elements missing."); return;
        }
        settingsModelInput.value = currentSettings.modelName;
        settingsTempSlider.value = currentSettings.temperature;
        settingsTempValueDisplay.textContent = currentSettings.temperature.toFixed(1);
        settingsMaxTokensInput.value = currentSettings.maxOutputTokens;
        settingsToolDropdown.value = currentSettings.enabledTool; // Set dropdown value
        settingsModalOverlay.classList.add('active');
    }
    function closeSettingsModal() {
        if (settingsModalOverlay) settingsModalOverlay.classList.remove('active');
    }
    function saveSettings() {
         if (!settingsModelInput || !settingsTempSlider || !settingsMaxTokensInput || !settingsToolDropdown) {
             console.error("Cannot save settings: Modal input elements missing."); return;
         }
         // Store old settings for comparison
         const oldSettings = {...currentSettings};

         // Read new values
         const newModelName = settingsModelInput.value.trim() || initialModelName;
         const newTemp = parseFloat(settingsTempSlider.value);
         const newMaxTokens = parseInt(settingsMaxTokensInput.value, 10);
         const newTool = settingsToolDropdown.value;

         const validatedMaxTokens = Math.max(2048, Math.min(newMaxTokens || 20000, 65536));

         // Update current settings
         currentSettings = {
             modelName: newModelName,
             temperature: newTemp,
             maxOutputTokens: validatedMaxTokens,
             enabledTool: newTool // Store selected tool
         };

         localStorage.setItem('chatSettings', JSON.stringify(currentSettings));
         updateSettingsDisplay();
         closeSettingsModal();
         console.log("Settings saved:", currentSettings);

         // Check if relevant settings changed to trigger new chat
         if (oldSettings.modelName !== currentSettings.modelName ||
             oldSettings.enabledTool !== currentSettings.enabledTool) {
             console.log("Model or Tool changed, starting new chat.");
             triggerNewChat(); // Call the new chat function
         } else if (oldSettings.temperature !== currentSettings.temperature ||
                    oldSettings.maxOutputTokens !== currentSettings.maxOutputTokens) {
             console.log("Temperature or Max Tokens changed. Settings applied on next message.");
             // Optionally notify user that these settings apply on next message in the *current* chat
             // Or force new chat for these too if desired: triggerNewChat();
         }
    }
    function loadSettings() {
         const savedSettings = localStorage.getItem('chatSettings');
         if (savedSettings) {
             try {
                 const parsedSettings = JSON.parse(savedSettings);
                 currentSettings.modelName = parsedSettings.modelName || initialModelName;
                 currentSettings.temperature = Math.max(0.0, Math.min(parseFloat(parsedSettings.temperature || 0.7), 2.0));
                 currentSettings.maxOutputTokens = Math.max(2048, Math.min(parseInt(parsedSettings.maxOutputTokens || 20000, 10), 65536));
                 currentSettings.enabledTool = parsedSettings.enabledTool || 'search'; // Load saved tool or default
             } catch (e) { console.error("Error parsing saved settings:", e); /* Use defaults below */ }
         }
         // Ensure defaults are set if nothing loaded or parsing failed
         currentSettings.modelName = currentSettings.modelName || initialModelName;
         currentSettings.temperature = currentSettings.temperature ?? 0.7;
         currentSettings.maxOutputTokens = currentSettings.maxOutputTokens || 20000;
         currentSettings.enabledTool = currentSettings.enabledTool || 'search';

         // Update modal inputs if they exist
         if(settingsModelInput) settingsModelInput.value = currentSettings.modelName;
         if(settingsTempSlider) settingsTempSlider.value = currentSettings.temperature;
         if(settingsTempValueDisplay) settingsTempValueDisplay.textContent = currentSettings.temperature.toFixed(1);
         if(settingsMaxTokensInput) settingsMaxTokensInput.value = currentSettings.maxOutputTokens;
         if(settingsToolDropdown) settingsToolDropdown.value = currentSettings.enabledTool; // Set dropdown default
         updateSettingsDisplay();
    }
    function updateSettingsDisplay() {
        if(headerModelDisplay) headerModelDisplay.textContent = currentSettings.modelName;
        if(footerModelDisplay) footerModelDisplay.textContent = currentSettings.modelName;
        if(footerTempDisplay) footerTempDisplay.textContent = `(Temp: ${currentSettings.temperature.toFixed(1)})`;
    }
    // Settings Event Listeners
    if (settingsButton) { settingsButton.addEventListener('click', openSettingsModal); }
    if (settingsCloseButton) { settingsCloseButton.addEventListener('click', closeSettingsModal); }
    if (settingsModalOverlay) { settingsModalOverlay.addEventListener('click', (event) => { if (event.target === settingsModalOverlay) closeSettingsModal(); }); }
    if (settingsSaveButton) { settingsSaveButton.addEventListener('click', saveSettings); }
    if (settingsTempSlider) { settingsTempSlider.addEventListener('input', () => { if(settingsTempValueDisplay) settingsTempValueDisplay.textContent = parseFloat(settingsTempSlider.value).toFixed(1); }); }

    // --- Password Modal Logic ---
    if (passwordForm) {
        passwordForm.addEventListener('submit', async (event) => {
            event.preventDefault();
            if (!passwordInput || !passwordSubmitButton || !passwordError) return;
            const enteredPassword = passwordInput.value;
            if (!enteredPassword) { passwordError.textContent = "Password cannot be empty."; return; }
            passwordError.textContent = "";
            passwordSubmitButton.disabled = true;
            passwordSubmitButton.textContent = "Checking...";
            try {
                const response = await fetch('/check_password', {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify({ password: enteredPassword })
                });
                const data = await response.json();
                if (response.ok && data.authenticated) {
                    isAuthenticated = true;
                    if (passwordModalOverlay) passwordModalOverlay.classList.remove('active');
                    document.body.classList.remove('password-protected');
                    disableChatInput(false);
                    if(promptInput) promptInput.focus();
                } else {
                    passwordError.textContent = data.error || "Incorrect password.";
                    passwordInput.value = "";
                    passwordInput.focus();
                }
            } catch (error) {
                console.error("Password check error:", error);
                passwordError.textContent = "Error checking password. Please try again.";
            } finally {
                if(passwordSubmitButton) {
                    passwordSubmitButton.disabled = false;
                    passwordSubmitButton.textContent = "Unlock";
                }
            }
        });
    } else if (window.passwordRequired) { console.error("Password form not found, but password protection is enabled!"); }

    // Function to disable/enable main chat input area
    function disableChatInput(disabled) {
        if (promptInput) promptInput.disabled = disabled;
        if (sendButton) sendButton.disabled = disabled;
        if (fileInput) fileInput.disabled = disabled;
        const footer = document.getElementById('main-footer');
        if (footer) {
            footer.style.opacity = disabled ? '0.5' : '1';
            footer.style.pointerEvents = disabled ? 'none' : 'auto';
        }
    }

    // --- New Chat Button ---
    if (newChatButton) { newChatButton.addEventListener('click', triggerNewChat); } // Use named function
    else { console.error("New chat button not found"); }

    // Function to handle starting a new chat (called by button or settings change)
    async function triggerNewChat() {
        console.log("Triggering New Chat...");
        try {
            const response = await fetch('/clear_chat', { method: 'POST' });
            if (!response.ok) { console.error("Failed to clear server-side chat session:", response.statusText); }
            else { console.log("Server-side chat session cleared."); }
        } catch (error) { console.error("Error calling /clear_chat endpoint:", error); }

        if (chatContainer) chatContainer.innerHTML = '';
        addInitialGreeting();
        clearFilePreviews();
        if (promptInput) { promptInput.value = ''; promptInput.style.height = 'auto'; }
        updateSendButtonState();
        console.log("Chat UI cleared.");
    }

    // --- Typing Indicator ---
    function showTypingIndicator() {
        if (typingIndicatorElement || !chatContainer) return;
        typingIndicatorElement = document.createElement('div');
        typingIndicatorElement.className = 'flex justify-start mb-4 typing-indicator-bubble';
        typingIndicatorElement.innerHTML = `
            <div class="max-w-xl lg:max-w-3xl p-3 rounded-lg shadow bg-model-bg-light dark:bg-model-bg-dark">
                <div class="flex items-center space-x-1 text-gray-600 dark:text-gray-400">
                    <span class="typing-dot"></span><span class="typing-dot"></span><span class="typing-dot"></span>
                </div>
            </div>`;
        chatContainer.appendChild(typingIndicatorElement);
        scrollToBottom();
    }
    function removeTypingIndicator() {
        if (typingIndicatorElement) { typingIndicatorElement.remove(); typingIndicatorElement = null; }
    }

    // --- Chat Submission ---
    if (chatForm) {
        chatForm.addEventListener('submit', async (event) => {
            event.preventDefault();
            if (!isAuthenticated && window.passwordRequired) { /* ... auth check ... */ return; }
            if (!promptInput || !sendButton || sendButton.disabled) return;
            const promptText = promptInput.value.trim();
            if (!promptText && attachedFiles.length === 0) return;

            // Pass copy of files for display
            displayMessage(promptText, [...attachedFiles], 'user');

            const formData = new FormData();
            formData.append('prompt', promptText);
            const filesToSend = [...attachedFiles];
            filesToSend.forEach(file => { formData.append('files', file, file.name); });
            // Send current settings
            formData.append('model_name', currentSettings.modelName);
            formData.append('temperature', currentSettings.temperature);
            formData.append('max_output_tokens', currentSettings.maxOutputTokens);
            formData.append('enabled_tool', currentSettings.enabledTool); // Send selected tool

            if(promptInput) { promptInput.value = ''; promptInput.style.height = 'auto'; }
            clearFilePreviews();

            updateSendButtonState();
            setLoading(true);
            showTypingIndicator();

            try {
                const response = await fetch('/chat', { method: 'POST', body: formData });
                removeTypingIndicator();
                setLoading(false);
                if (!response.ok) {
                    const errorData = await response.json().catch(() => ({ error: "Failed to parse error response." }));
                    console.error('Error response from server:', errorData);
                    displayMessage(`Server Error: ${errorData.error || response.statusText}`, [], 'error', null, []);
                    return;
                }
                const data = await response.json();
                // Update displayed model name if backend used a different one (e.g., fallback)
                if (data.model_used && headerModelDisplay && footerModelDisplay) {
                     headerModelDisplay.textContent = data.model_used;
                     footerModelDisplay.textContent = data.model_used;
                     // Optionally update currentSettings.modelName if you want the UI change to persist
                     // currentSettings.modelName = data.model_used;
                     // localStorage.setItem('chatSettings', JSON.stringify(currentSettings));
                }
                displayMessage(null, [], 'model', data.finish_reason, data.parts);
            } catch (fetchError) {
                console.error('Fetch error:', fetchError);
                removeTypingIndicator();
                setLoading(false);
                displayMessage(`Network error: ${fetchError.message}`, [], 'error', null, []);
            }
        });
    } else { console.error("Chat form not found"); }

    // --- Display Messages ---
    function addInitialGreeting() {
         if (!chatContainer) return;
         if (chatContainer.children.length === 0 && isAuthenticated) {
              const greetingText = "Hello! How can I help you today?";
              const messageDiv = document.createElement('div');
              messageDiv.className = `flex justify-start mb-4`;
              const bubbleDiv = document.createElement('div');
              bubbleDiv.className = `message-bubble max-w-xl lg:max-w-3xl p-3 rounded-lg shadow bg-model-bg-light dark:bg-model-bg-dark`;
              bubbleDiv.innerHTML = `<p>${greetingText}</p>`;
              messageDiv.appendChild(bubbleDiv);
              chatContainer.appendChild(messageDiv);
         }
    }

    function displayMessage(text, files = [], role, finishReason = null, parts = []) {
         if (!chatContainer) { console.error("Cannot display message: chatContainer element not found."); return; }

         const messageDiv = document.createElement('div');
         messageDiv.className = `flex ${role === 'user' ? 'justify-end' : 'justify-start'} mb-4`;
         const bubbleDiv = document.createElement('div');
         bubbleDiv.className = `message-bubble max-w-xl lg:max-w-3xl p-3 rounded-lg shadow overflow-hidden ${role === 'user' ? 'bg-user-bg-light dark:bg-user-bg-dark' : 'bg-model-bg-light dark:bg-model-bg-dark'}`;
         let contentHTML = '';

         // --- Render based on role ---
         if (role === 'user') {
             // User messages use 'text' and 'files' arguments
             if (files.length > 0) {
                 contentHTML += '<div class="flex flex-wrap gap-2 mb-2">';
                 files.forEach(file => {
                      let iconClass = 'fa-file';
                     if (file.type.startsWith('image/')) iconClass = 'fa-file-image';
                     else if (file.type === 'application/pdf') iconClass = 'fa-file-pdf';
                     else if (file.type.startsWith('text/')) iconClass = 'fa-file-alt';
                     else if (file.type.includes('python') || file.name.endsWith('.py')) iconClass = 'fa-file-code';
                     contentHTML += `<div class="bg-gray-200 dark:bg-gray-600 p-1 px-2 rounded text-xs flex items-center gap-1"> <i class="fas ${iconClass} text-gray-600 dark:text-gray-400"></i> <span>${file.name}</span> </div>`;
                 });
                 contentHTML += '</div>';
             }
             if (text) {
                 const escapedText = text.replace(/</g, "&lt;").replace(/>/g, "&gt;");
                 contentHTML += `<p class="whitespace-pre-wrap">${escapedText}</p>`;
             } else if (files.length > 0 && !text) {
                 contentHTML += `<p class="italic text-sm text-gray-500 dark:text-gray-400">[Uploaded ${files.length} file(s)]</p>`;
             }
         } else { // Model or Error roles - process the 'parts' array
             if (parts && parts.length > 0) {
                 // console.log("Rendering received parts:", parts); // Optional log
                 parts.forEach(part => {
                     // console.log("Processing part:", part); // Optional log
                     switch (part.type) {
                         case 'text':
                             if (typeof marked !== 'undefined' && typeof DOMPurify !== 'undefined') {
                                 const unsafeHTML = marked.parse(part.content || ""); // Use configured marked
                                 const safeHTML = DOMPurify.sanitize(unsafeHTML);
                                 contentHTML += `<div class="prose dark:prose-invert max-w-none">${safeHTML}</div>`;
                             } else {
                                 contentHTML += `<p class="whitespace-pre-wrap">${(part.content || "").replace(/</g, "&lt;").replace(/>/g, "&gt;")}</p>`; // Fallback
                             }
                             break;
                         case 'executable_code':
                             const langClass = part.language ? `language-${part.language}` : 'language-plaintext';
                             const escapedCode = (part.code || "").replace(/</g, "&lt;").replace(/>/g, "&gt;");
                             // Wrap manually for structure before highlighting
                             contentHTML += `<div class="code-container my-2"><pre><code class="${langClass}">${escapedCode}</code></pre></div>`;
                             break;
                         case 'code_result':
                             const escapedOutput = (part.output || "").replace(/</g, "&lt;").replace(/>/g, "&gt;");
                             contentHTML += `<pre class="code-output bg-gray-100 dark:bg-gray-900 p-2 rounded text-sm whitespace-pre-wrap my-2"><strong>[Result: ${part.outcome || 'UNKNOWN'}]</strong>\n${escapedOutput}</pre>`;
                             break;
                         case 'image':
                             // console.log("Rendering image part with URL:", part.url); // Optional log
                             contentHTML += `<img src="${part.url}" alt="Generated Image" class="max-w-full h-auto rounded my-2 block" onerror="this.alt='Error loading image'; this.style.display='none'; console.error('Error loading image:', this.src);">`; // Added onerror
                             break;
                          case 'search_suggestions': // Handle search suggestions
                             // IMPORTANT: Insert HTML directly without sanitization
                             contentHTML += `<div class="search-suggestions my-3">${part.html_content || ''}</div>`;
                             break;
                         default:
                             console.warn("Unknown part type received:", part.type);
                             const escapedContent = JSON.stringify(part.content || "").replace(/</g, "&lt;").replace(/>/g, "&gt;");
                             contentHTML += `<p class="text-xs text-red-500">[Unsupported content type: ${part.type}]</p><pre class="text-xs">${escapedContent}</pre>`;
                     }
                 });
             } else if (text) { // Fallback for simple error messages passed as text
                 contentHTML += `<p class="whitespace-pre-wrap">${text.replace(/</g, "&lt;").replace(/>/g, "&gt;")}</p>`;
             } else {
                 contentHTML += `<p>[No content received]</p>`;
             }

             // Finish Reason (Append after processing parts)
             if (finishReason && finishReason !== 'STOP' && finishReason !== 'MAX_TOKENS' && finishReason !== 'FINISH_REASON_UNSPECIFIED') {
                 contentHTML += `<p class="text-xs text-gray-500 dark:text-gray-400 mt-2 border-t border-border-light dark:border-border-dark pt-1">Finish Reason: ${finishReason}</p>`;
             }
         }

         // Error Styling
         if (role === 'error') {
              bubbleDiv.classList.add('bg-red-100', 'dark:bg-red-900', 'text-red-700', 'dark:text-red-300');
         }

         // Set bubble content
         bubbleDiv.innerHTML = contentHTML;
         messageDiv.appendChild(bubbleDiv);
         chatContainer.appendChild(messageDiv);

         // Apply Syntax Highlighting and Add Copy Buttons AFTER rendering
         if (role === 'model' || role === 'error') {
            try {
                if (typeof hljs !== 'undefined') {
                    highlightCodeInElement(bubbleDiv);
                }
                addCopyButtonsToCodeBlocks(bubbleDiv);
            } catch(e) {
                console.error("Error applying highlighting or copy buttons:", e);
            }
         }

         scrollToBottom();
    }

    // --- Apply Syntax Highlighting ---
    function highlightCodeInElement(containerElement) {
        if (typeof hljs === 'undefined') { return; } // Guard clause
        // Find code blocks that haven't been highlighted yet
        const codeBlocks = containerElement.querySelectorAll('pre code:not(.hljs)');
        codeBlocks.forEach((codeElement) => {
            try {
                 hljs.highlightElement(codeElement);
            } catch (error) {
                console.error("Error highlighting element:", error, codeElement);
            }
        });
    }


    // --- Add Copy Buttons ---
    function addCopyButtonsToCodeBlocks(containerElement) {
        // Add check for clipboard API support
        if (typeof navigator.clipboard === 'undefined' || typeof navigator.clipboard.writeText === 'undefined') {
             console.warn("Clipboard API (writeText) not available, copy buttons will not function correctly.");
             // return; // Decide whether to add non-functional buttons or not
        }
        // Find PRE elements that likely contain code (have a CODE child with a language class)
        const codeBlocks = containerElement.querySelectorAll('pre code[class*="language-"]');
        codeBlocks.forEach(codeElement => {
            const preElement = codeElement.parentNode; // Get the parent <pre>
            if (!preElement || preElement.tagName !== 'PRE') return; // Ensure parent is <pre>
            // Check if already wrapped
            if (preElement.parentNode.classList.contains('code-block-wrapper')) return;

            const wrapper = document.createElement('div');
            wrapper.className = 'code-block-wrapper'; // Class for relative positioning

            const button = document.createElement('button');
            const copyIconHtml = '<i class="far fa-copy mr-1"></i>'; // Store icon HTML
            button.innerHTML = `${copyIconHtml}Copy`;
            // Tailwind classes for styling
            button.className = 'copy-code-button bg-gray-300 dark:bg-gray-700 text-gray-800 dark:text-gray-200 hover:bg-gray-400 dark:hover:bg-gray-600';
            button.style.fontFamily = 'inherit'; // Match surrounding font
            button.setAttribute('title', 'Copy code'); // Add tooltip

            button.addEventListener('click', (e) => {
                e.stopPropagation(); // Prevent potential parent listeners
                const codeToCopy = codeElement.textContent || "";

                // Check clipboard API availability *inside* the handler as well
                if (!navigator.clipboard || !navigator.clipboard.writeText) {
                    console.error("Clipboard API not available or not permitted.");
                    button.innerHTML = '<i class="fas fa-times mr-1"></i> Error';
                     button.classList.add('bg-red-500', 'text-white');
                     button.classList.remove('bg-gray-300', 'dark:bg-gray-700', 'text-gray-800', 'dark:text-gray-200');
                     button.disabled = true;
                     setTimeout(() => { // Reset after timeout
                         button.innerHTML = `${copyIconHtml}Copy`;
                         button.classList.remove('bg-red-500', 'text-white');
                         button.classList.add('bg-gray-300', 'dark:bg-gray-700', 'text-gray-800', 'dark:text-gray-200');
                         button.disabled = false;
                     }, 3000);
                    return; // Stop execution
                }

                // Try to copy
                navigator.clipboard.writeText(codeToCopy).then(() => {
                    button.innerHTML = '<i class="fas fa-check mr-1"></i>Copied!'; // Checkmark icon
                    // Success state styling
                    button.classList.add('bg-green-500', 'dark:bg-green-600', 'text-white');
                    button.classList.remove('bg-gray-300', 'dark:bg-gray-700', 'text-gray-800', 'dark:text-gray-200', 'hover:bg-gray-400', 'dark:hover:bg-gray-600');
                    button.disabled = true;

                    setTimeout(() => {
                        button.innerHTML = `${copyIconHtml}Copy`; // Restore icon + text
                        // Restore original styling
                        button.classList.remove('bg-green-500', 'dark:bg-green-600', 'text-white');
                        button.classList.add('bg-gray-300', 'dark:bg-gray-700', 'text-gray-800', 'dark:text-gray-200');
                         button.disabled = false;
                    }, 2000);
                }).catch(err => {
                    console.error('Failed to copy code via navigator.clipboard:', err); // Log specific error
                    button.innerHTML = '<i class="fas fa-times mr-1"></i>Failed!'; // More specific error text
                    // Error styling
                    button.classList.add('bg-red-500', 'text-white');
                    button.classList.remove('bg-gray-300', 'dark:bg-gray-700', 'text-gray-800', 'dark:text-gray-200');
                    button.disabled = false; // Allow retry on error

                    setTimeout(() => {
                         button.innerHTML = `${copyIconHtml}Copy`; // Restore icon + text
                         button.classList.remove('bg-red-500', 'text-white');
                         button.classList.add('bg-gray-300', 'dark:bg-gray-700', 'text-gray-800', 'dark:text-gray-200');
                    }, 3000); // Longer timeout for error message
                });
            });

            // Wrap the <pre> and add the button
            if (preElement.parentNode) {
                preElement.parentNode.insertBefore(wrapper, preElement);
                wrapper.appendChild(preElement); // Move <pre> inside wrapper
                wrapper.appendChild(button);    // Add button inside wrapper
            } else { console.warn("Could not wrap pre element, parentNode missing."); }
        });
    }


     // --- Utility Functions ---
    function setLoading(isLoading) {
         if (!sendButton || !sendIcon || !loadingIndicator) return;
         if (isLoading) {
             sendIcon.classList.add('hidden');
             loadingIndicator.classList.remove('hidden');
             sendButton.disabled = true;
             sendButton.classList.add('cursor-wait');
         } else {
             sendIcon.classList.remove('hidden');
             loadingIndicator.classList.add('hidden');
             sendButton.classList.remove('cursor-wait');
             updateSendButtonState(); // Re-evaluate based on content
         }
     }
    function scrollToBottom() {
         if(chatContainer) chatContainer.scrollTop = chatContainer.scrollHeight;
     }

    // --- Enter Key Listener ---
    if (promptInput) {
        promptInput.addEventListener('keydown', (event) => {
            if (event.key === 'Enter' && !event.shiftKey) {
                event.preventDefault(); // Prevent newline
                // Trigger button click only if it's not disabled
                if (sendButton && !sendButton.disabled) {
                    sendButton.click(); // Programmatically click the button
                }
            }
        });
    } else { console.error("Prompt input not found for keydown listener"); }

}); // End DOMContentLoaded
